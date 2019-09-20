function [counts, covarout] = covar_annot(Mu, covars, T)

%
%process mutation file
demand_fields(Mu, {'chr' 'pos' 'count' 'context'});

Mu = sort_struct(Mu, {'chr' 'pos'});
[cu, cui] = unique(Mu.chr);

contexts = unique(Mu.context);

%
%handle case of only processing a single covariate
if ischar(covars), covars = {covars}; end

counts = cell(length(cu), 1);
covarout = cell(length(cu), length(covars));

%
%load covariates into memory
C = cell(length(covars), 1);
for c = 1:length(covars),
  C{c} = [];
  C{c}.file = direc(covars{c});
  C{c}.chr = convert_chr(regexprep(C{c}.file, '.*/chr([0-9XY]+)\.mat', '$1'));
  C{c} = reorder_struct(C{c}, ~isnan(C{c}.chr));

  if slength(C{c}) ~= 24, error(['Invalid number of chromosomes in covariate directory "' covars{c} '"!']); end

  C{c} = sort_struct(C{c}, 'chr');
end

%
%process restriction interval list
if exist('T', 'var'),
  if ~isstruct(T), error('Restriction interval list must be specified as struct!'); end
  demand_fields(T, {'chr' 'start' 'end'})
  T = sort_struct(T, {'chr' 'start' 'end'});

  [~, rcui] = unique(T.chr);
  tbdy = [rcui [rcui(2:end) - 1; slength(T)]];
end

%loop over chromosomes
for x = [cui [cui(2:end) - 1; slength(Mu)] cu (1:length(cu))']',
  i = x(1); j = x(2); k = x(3); l = x(4);

  %loop over covariates
  for c = 1:length(covars), 
    Q = load(C{c}.file{k});
    demand_fields(Q, {'covu' 'covar_mat'});

    %subset covariates to contexts
    %NB: we use this awkward trick to sum across columns because Matlab's sum() function 
    %    reconstitutes sparse matrices as full, and then re-sparsifies them -- SLOW!
    tmp = Q.covar_mat(:, contexts);
    covar_mat_subset = tmp(:, 1);
    for q = 2:size(tmp, 2),
      covar_mat_subset = covar_mat_subset + tmp(:, q);
    end
    clear tmp

    %restrict territory (if specified)
    if exist('T', 'var'),
      bdy = tbdy(k, :);
      is_ok = false(size(covar_mat_subset));

      for y = bdy(1):bdy(2),
	is_ok(T.start(y):T.end(y)) = true;
      end

      covar_mat_subset = covar_mat_subset.*is_ok;
    end

    %index covariates for mutated positions
    midx = full(covar_mat_subset(Mu.pos(i:j)));

    if any(midx == 0),
      error(['Covariate domain mismatch -- are you trying to annotate noncoding ' ...
	     'mutations with coding covariates, or non-restricted mutations with a ' ...
             'restriction interval list?']);
    end

    %index covariates for non-mutated positions
    zidx = covar_mat_subset > 0;
    zidx(Mu.pos(i:j)) = 0;
    zidx = full(covar_mat_subset(zidx));

    %concatenate covariates for mutated positions with covariates for non-mutated positions
    nm = length(midx); nz = length(zidx);
    covarout{l, c} = NaN(nm + nz, 1);
    covarout{l, c}(1:nm) = Q.covu(midx);
    covarout{l, c}((nm + 1):end) = Q.covu(zidx);

    %output counts with same concatenation
    if c == 1, counts{l} = [Mu.count(i:j); zeros(nz, 1)]; end

    fprintf('(%d/%d, %d/%d) ', c, length(covars), l, length(cu));
  end
end

%concatenate 
counts = cat(1, counts{:});
covarout = cell2mat(covarout);
