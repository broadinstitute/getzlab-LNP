function covar_index(context_dir, targlist, covar_file, outdir, ztrans)

usage = 'covar_index(<directory containing contexts for each chromosome>, <target list>, <covariates>, <output directory>, <Z transform flag>)';

%
%argument handling
if nargin == 4,
  ztrans = 0;
elseif nargin == 5 && (islogical(ztrans) || ztrans == 1 || ztrans == 0),
else
  error(['Usage: ' usage]);
end

%
%ensure directory containing contexts is OK
C = [];
C.file = direc(context_dir);
C.chr = convert_chr(regexprep(C.file, '.*/chr([0-9XY]+)\.mat', '$1'));
C = reorder_struct(C, ~isnan(C.chr));

if slength(C) ~= 24, error('Invalid number of chromosomes in context track directory!'); end

C = sort_struct(C, 'chr');

%
%ensure output directory is OK
ede(outdir)

%
%load and process target list

%TODO: check that fields are numeric; accept as both struct or TSV
T = targlist; clear targlist
demand_fields(T, {'chr' 'start' 'end'});

T = sort_struct(T, {'chr' 'start' 'end'});
[cu, cui] = unique(T.chr);

%
%load and process covariates list

%we infer whether the covariates are positionwise encoded (i.e. FWB format) or intervalwise encoded
%(i.e. MAT format) based on the file extension.

%MAT input must be formatted like a target list
%TODO: also allow TSV or struct input for this

if ischar(covar_file),
  if grepmi('\.fwb$', {covar_file}),
    Q = org.broadinstitute.cga.tools.seq.FixedWidthBinary(covar_file);
    Q.setNullVal(0); 
  elseif grepmi('\.mat$', {covar_file}),
    tmp = load(covar_file);
    if ~isfield(tmp, 'Q'), error('Variable in MAT file must be named "Q"'); end
    Q = tmp.Q; clear tmp 
  else
    error('Covariate input file must either be FWB or MAT formatted!');
  end
end

if isjava(Q), %nothing to do, since we presume that FWB class was instantiated in our code 
              %XXX: this will not work if a Java class was provided directly
elseif isstruct(Q),
  demand_fields(Q, {'chr' 'start' 'end' 'covar'});
  Q = sort_struct(Q, {'chr' 'start' 'end'});

  %TODO: check for overlapping intervals
else
  error('Invalid format for covariates!')
end

%
%variables for storing incremental mean/variance for Z-standardization
count_z = 0;
weight_sum = NaN;
covar_mean = NaN;
covar_std = NaN;

covus = cell(length(cui), 1);
covar_mats = cell(length(cui), 1);

%
%loop over chromosomes
for x = [cui [cui(2:end) - 1; slength(T)] cu]'
  i = x(1); j = x(2); k = x(3);

  %load pentamer track for this chromosome
  tmp = load(C.file{k});
  if ~isfield(tmp, 'categ'), error('Context track variable must be named "categ"!'); end
  categ = tmp.categ; clear tmp

  %mask targets
  is_coding = false(size(categ)); 
  for l = i:j,
    is_coding(T.start(l):T.end(l)) = true;
  end
  codposes = find(is_coding);

  %load covariate values, either from FWB track or from interval list
  ctrack = NaN(nnz(is_coding), 1);
  if isjava(Q),
    ctrack = double(Q.get(k, codposes));
  else
    targ = map_mutations_to_targets_fast([k*ones(size(codposes)) codposes], [Q.chr Q.start Q.end]);

    idx = targ > 0;
    ctrack(idx) = Q.covar(targ(idx));
    ctrack(~idx) = -Inf; %standin for NaN, since NaN doesn't get properly uniqued
  end
  ctrack(isnan(ctrack)) = Inf; %-Inf == not in covariate track; +Inf == NaN in covariate track

  %Z-standardizing covariates (if requested) requires iteratively calculating mean/variance
  if ztrans && ~all(isinf(ctrack)),
    idx = ~isinf(ctrack);
    n = nnz(idx);
    chrmean = mean(ctrack(idx));
    chrvar = var(ctrack(idx), 1);

    if count_z == 0,
      weight_sum = n;
      covar_mean = chrmean;
      covar_std = chrvar + chrmean^2;
    else
      weight_sum = weight_sum + n;
      covar_mean = covar_mean + n/weight_sum*(chrmean - covar_mean);
      covar_std = covar_std + n/weight_sum*(chrvar + chrmean^2 - covar_std);
    end

    count_z = 1;
  end

  %generate covariate LuT and contextwise pointers
  [covus{k}, ~, covuj] = unique(ctrack); 
  covar_mats{k} = sparse(codposes, categ(is_coding), covuj, size(categ, 1), 1024);

  fprintf('%d ', k);
end

for k = 1:length(cui), 
  if ztrans,
    covu = (covus{k} - covar_mean)./sqrt(covar_std - covar_mean^2);
  else
    covu = covus{k};
  end
  covar_mat = covar_mats{k};

  save([outdir '/chr' num2str(k) '.mat'], 'covu', 'covar_mat'); 
end
