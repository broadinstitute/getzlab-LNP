function pois_LN_reg_wrapper(M, ch, chfield, basechange, P)

if isdeployed,
  ch = str2double(strsplit(ch, ','));
  basechange = str2double(basechange);
end

%handle parameters
if exist('P', 'var') && ~isempty(P) && ischar(P),
  if all(P((end - 2):end) == 'txt'),
    P = process_params_file([], P);
  elseif all(P((end - 2):end) == 'mat'),
    tmp = load(P);
    if ~isfield(tmp, 'P'), error('Parameters struct must be named "P"!'); end
    if ~isstruct(tmp.P), error('Parameters must be specified as a struct!'); end
    P = tmp.P; clear tmp
  else
    error('Unknown format for params file!');
  end
elseif exist('P', 'var') && isstruct(P),
else
  P = [];
end

P = impose_default_value(P, 'outdir', '.');
P = impose_default_value(P, 'categs_file', 'ref/context_1025_categs.txt');
P = impose_default_value(P, 'covar_dir', 'ref/cod_mutsig_covars');
P = impose_default_value(P, 'terr_file', 'ref/cod_c512_covars.align75_filtered/v1');
P = impose_default_value(P, 'interval_list', 'none');
P = impose_default_value(P, 'categ_regex', '(.) in (..)_(..)');
P = impose_default_value(P, 'categ_sub', '$2$1$3');
P = impose_default_value(P, 'pool_basechanges', true);
P = impose_default_value(P, 'tailor_hyperparameters', 0);
P = impose_default_value(P, 'init_from_MAP', 0);
P = impose_default_value(P, 'MAP_NR_iterations', 20);
P = impose_default_value(P, 'min_MAP_gradient_norm', 1e-8);
P = impose_default_value(P, 'MAP_only', 0);
P = impose_default_value(P, 'compute_posterior_predictive', 0);
P = impose_default_value(P, 'pp_recurrence_floor', 0);
P = impose_default_value(P, 'fnum', bitor(ch(1), 1024));
P = impose_default_value(P, 'tier', 4);
P = impose_default_value(P, 'annotate_and_exit', false);
P = impose_default_value(P, 'skip_figures', true);
P = impose_default_value(P, 'save_output', 1);
P = impose_default_value(P, 'data_downsamp', Inf);

%can either be scalar (representing overall cohort size), or vector (for weighting each site by territory)
P = impose_default_value(P, 'log_exposure', 0);

%load and subset mutations
if ischar(M),
  demand_file(M)

  if exist(M, 'dir') == 7, %assume it's a saved M struct
    M = loadM(M); %TODO: need more safety checks here.
    M.mut.gene = M.gene.name(M.mut.gene_idx);
    Mu = M.mut; clear M
  elseif exist(M, 'file') == 2, %assume it's a TSV file; convert appropriately
    Mu = load_struct(M);

    % make sure required fields exist (renaming if necessary)
    rf = {
	{'chr', 'Chromosome'},
	{'pos', 'Position', 'start', 'Start_position'}
    };
  end
elseif isstruct(M),
  Mu = M; clear M
else
  error('Invalid format for M!');
end
demand_fields(Mu, {'chr' 'pos' chfield 'is_coding' 'tier'});

%XXX: for now, we assume no multilevel modeling -- we either operate on sums or independent changes
if P.pool_basechanges, %i.e. use the sum
  demand_field(Mu, 'count');
else
  demand_field(Mu, 'count_nb');
  Mu.count = Mu.count_nb(:, basechange);
end

Mu = reorder_struct(Mu, ismember(Mu.tier, P.tier) & Mu.is_coding & ismember(Mu.(chfield), ch));
Mu = sort_struct(Mu, {'chr' 'pos'});

%map mutations to targets if interval list is given
if ~strcmp(P.interval_list, 'none'),
  T = load_struct(P.interval_list);
  demand_fields(T, {'chr' 'start' 'end'})
  T = makeapn(T);
  T = sort_struct(T, {'chr' 'start' 'end'});

  targ = map_mutations_to_targets_fast([Mu.chr Mu.pos], [T.chr T.start T.end]);
  Mu = reorder_struct(Mu, targ > 0);
end

Mu = rename_field(Mu, chfield, 'context');

Mu.p = NaN(slength(Mu), 1);

%load context names
categs = load_struct(P.categs_file);
categs.name = regexprep(categs.name, P.categ_regex, P.categ_sub);

%if we're dealing with a single base change, then add the base change to categs
if ~P.pool_basechanges,
  B = 'ACGT';
  b('ACGT') = 1:4;
  nbidx = [2 3 4; 1 3 4; 1 2 4; 1 2 3];

  x = cat(1, categs.name{ch}); x = b(x(:, 3));
  bc = basechange*ones(size(ch)); bc(x > 2) = 4 - bc(x > 2);
  nbs = nbidx(sub2ind(size(nbidx), x, bc));

  for x = [ch; nbs],
    i = x(1); j = x(2);
 
    categs.name{i} = [categs.name{i}(1:2) '(' categs.name{i}(3) '->' B(j) ')' categs.name{i}(4:5)];
  end
end

%set default output name based on categ names
P = impose_default_value(P, 'output_name', [strjoin(num2cellstr(ch), ',') '-' strjoin(categs.name(ch), ',')]);

%load covariates
if ~strcmp(P.covar_dir, 'none'),
  if ~strcmp(P.interval_list, 'none'),
    [M C] = covar_annot(Mu, P.covar_dir, T);
  else
    [M C] = covar_annot(Mu, P.covar_dir);
  end 
else %no covariates
  if ~strcmp(P.interval_list, 'none'),
    M = covar_annot(Mu, P.terr_file, T);
  else
    M = covar_annot(Mu, P.terr_file);
  end
  C = zeros(size(M));
end
Mu = reorder_struct(Mu, Mu.count > 0);

%downsample data (if requested)
if ~isinf(P.data_downsamp),
  %P.data_downsamp can either be expressed as a number or factor
  %if the latter, convert to a number
  if P.data_downsamp > 0 && P.data_downsamp < 1,
    P.data_downsamp = round(P.data_downsamp*size(M, 1));
  end

  %binary index of sites to sample
  ds_idx = false(size(M));
  ds_idx(randsample(size(M, 1), P.data_downsamp)) = true;

  %index M -> Mu, since we will also downsample Mu struct accordingly
  Mu_idx = zeros(size(M));
  Mu_idx(find(M > 0)) = 1:nnz(M > 0);

  %downsample M/C
  M = M(ds_idx);
  C = C(ds_idx, :);

  %downsample Mu
  Mu = reorder_struct(Mu, Mu_idx(ds_idx & Mu_idx > 0));
end

%handle categorical covariates
%we do this after downsampling to ensure that indices are contiguous
if ~strcmp(P.covar_dir, 'none'),
  %handle fixed effect categorical covariates -- convert numerical categories to binary matrix
  if isfield(P, 'fixed_effect_cat_covar_index'),
    %TODO: add error checking to ensure this index is OK
    %TODO: allow for multiple binary fields.  for now, we only allow one
    if strcmp(P.fixed_effect_cat_covar_index, 'last'), P.fixed_effect_cat_covar_index = size(C, 2); end 
    if length(P.fixed_effect_cat_covar_index) > 1 || P.fixed_effect_cat_covar_index > size(C, 2),
      error('At the moment, you can only specify a single covariate to be categorical.')
    end

    %remove categories that have no mutations
    h = sparse(C(:, P.fixed_effect_cat_covar_index), M + 1, 1);
    zidx = ismember(C(:, P.fixed_effect_cat_covar_index), find(sum(h(:, 2:end), 2) == 0));
    M(zidx) = [];
    C(zidx, :) = [];

    %generate dummy variablesq
    [cov_u, ~, cov_uj] = unique(C(:, P.fixed_effect_cat_covar_index));
    
    C(:, P.fixed_effect_cat_covar_index) = [];
    C = [C full(sparse(1:size(C, 1), cov_uj, 1))];

    %avoid the dummy variable trap!
    %need to remove one of the binary covariates or else the model is overparameterized.
    C(:, end) = [];
    cov_u(end) = [];

    %ensure the names still sync up
    if isfield(P, 'covar_names'),
      P.covar_names(P.fixed_effect_cat_covar_index) = [];
      P.covar_names = [P.covar_names; num2cellstr(cov_u)];
    else
      P.covar_names = num2cellstr(cov_u);
    end
  end

  %handle random effect categorical covariates -- ensure they're uniqued and zero indexed
  if isfield(P, 'random_effect_cat_covar_index'),
    %TODO: check if this conflicts with P.fixed_effect_cat_covar_index
    if strcmp(P.random_effect_cat_covar_index, 'last'),
      P.random_effect_cat_covar_index = size(C, 2);
    end 
    if length(P.random_effect_cat_covar_index) > 1 || P.random_effect_cat_covar_index > size(C, 2),
      error('At the moment, you can only specify a single covariate to be categorical.')
    end

    %remove categories that have no mutations
    h = sparse(C(:, P.random_effect_cat_covar_index), M + 1, 1);
    zidx = ismember(C(:, P.random_effect_cat_covar_index), find(sum(h(:, 2:end), 2) == 0));
    M(zidx) = [];
    C(zidx, :) = [];

    %make indicators contiguous
    [~, ~, cov_uj] = unique(C(:, P.random_effect_cat_covar_index));
    C(:, P.random_effect_cat_covar_index) = cov_uj;
  end
end

fprintf('\n');

%calculations start here
fprintf('\nch: %d\n\n', ch);

%scaled pois. MLE used as lambda hyperparameter
P.sepsi = 3.5;
P.use_tailored_MH = true;

%use glmfit for initial parameter estimates
%whether we include an intercept term depends on whether we have categorical random effects.
C(isnan(C)) = 0; C(isinf(C)) = 0;
if isfield(P, 'random_effect_cat_covar_index'),
  const = 'off';
  C2 = C;
  C2(:, P.random_effect_cat_covar_index) = [];

  C2 = [C2 full(sparse(1:size(C2, 1), ...
                       C(:, P.random_effect_cat_covar_index), 1))];
else
  const = 'on';
  C2 = C;
end
BetaP = glmfit(C2, M, 'poisson', 'offset', ones(size(M)).*P.log_exposure, 'constant', const);

%tailor hyperparameters (if requested)
if P.tailor_hyperparameters,
  if isfield(P, 'random_effect_cat_covar_index'),
    P.mumu = BetaP(size(C, 2):end);
    P.mubeta = BetaP(1:(size(C, 2) - 1));
  else
    P.mumu = BetaP(1);
    P.mubeta = BetaP(2:end);
  end

  %TODO: some heuristic for tau hyperparameters?
end

%calculate optimal values of mu/tau/beta
if P.init_from_MAP,
  PP = P;
  PP.niter = P.MAP_NR_iterations;
  PP.min_gradient_norm = P.min_MAP_gradient_norm;
  PP.m0 = BetaP(1);
  PP.beta0 = BetaP(2:end);

  %hyperparameters
  PP = impose_default_value(PP, 'taua', 1);
  PP = impose_default_value(PP, 'taub', 10);
  PP = impose_default_value(PP, 'mumu', 0);
  PP = impose_default_value(PP, 'mutau', 1);
  PP = impose_default_value(PP, 'invcovmu', 0.5*eye(size(C, 2)));
  PP = impose_default_value(PP, 'mubeta', zeros(size(C, 2), 1));

  [m0 t0 beta0] = regr_NR(M, P.log_exposure, C, PP, PP);

  P.m0 = m0(end);
  P.t0 = t0(end);
  P.beta0 = beta0(end, :)';

  P.post_abs_mu = P.m0;
  P.post_abs_tau = P.t0;
  P.post_abs_Beta = P.beta0;
else
  P = impose_default_value(P, 't0', 1);

  if isfield(P, 'random_effect_cat_covar_index'),
    P = impose_default_value(P, 'm0', BetaP(size(C, 2):end));
    P = impose_default_value(P, 'beta0', BetaP(1:(size(C, 2) - 1)));
  else
    P = impose_default_value(P, 'm0', BetaP(1));
    P = impose_default_value(P, 'beta0', BetaP(2:end));
  end

  %if we are running sans covariates, set beta0 to 0
  if size(P.beta0, 1) == 0, P.beta0 = 0; end
end

if P.annotate_and_exit,
  save([P.outdir '/' P.output_name '_annot-only.mat'], 'Mu', 'M', 'C', 'P')
  return;
end

if P.MAP_only,
  if ~P.init_from_MAP, error('Cannot output MAP if we did not calculate it!'); end

  if ~P.save_output, error('Nothing to do!'); return; end

  save([P.outdir '/' P.output_name '_MAPonly.mat'], 'beta0', 'm0', 't0', 'Mu', 'M', 'C', 'P');
  return;
end

%run regression
[Beta, epsi, mu, tau, logmarg_lik] = pois_LN_reg(M, P.log_exposure, C, P);
post_pred = NaN;
pvalues = NaN;
pp_idx = NaN;

%compute posterior predictive
if P.compute_posterior_predictive,
  if P.pp_recurrence_floor > 0, %compute for loci with recurrence above floor
    idx = M >= P.pp_recurrence_floor;

    [post_pred pvalues] = regr_post_pred(mu, tau, Beta, M(idx), P.log_exposure, C(idx, :), P);

    %annotate MAF with mean p_max
    Mu.p(Mu.count >= P.pp_recurrence_floor) = pvalues(:, 3, 1);
  else %compute for all loci
    [post_pred pvalues] = regr_post_pred(mu, tau, Beta, M, P.log_exposure, C, P);

    %annotate MAF with mean p_max
    Mu.p = pvalues(M > 0, 3, 1);
  end 
else
  post_pred = NaN;
  pvalues = NaN;
end

%save output
if P.save_output,
  save([P.outdir '/' P.output_name '.mat'], 'Beta', 'epsi', 'mu', 'tau', 'logmarg_lik', 'logmarg_lik_unif', 'post_pred', 'pvalues', 'Mu', 'M', 'C', 'P')
end

if P.skip_figures,
  return;
end

fnum = P.fnum;
regr_figure(M, Mu, C, post_pred, pvalues, mu, tau, Beta, fnum, P);

if P.save_output,
  set(fnum, 'PaperPositionMode', 'auto')
  print([P.outdir '/' P.output_name '.png'], '-dpng')

  if isdeployed,
    close(fnum)
  end
end
