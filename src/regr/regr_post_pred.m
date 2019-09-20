function [post_pred pvalues] = regr_post_pred(mu, tau, Beta, xpp, ypp, cpp, P)

if ~exist('P', 'var'),
  error('Parameters must be provided in P struct!');
end

demand_fields(P, {'burnin' 'niter'})

P = impose_default_value(P, 'pp_buffer_size', 100);
P = impose_default_value(P, 'thin', 3);

%look up Hermite polynomial roots/quadrature weights
X = [];
[X.h_r X.q_w] = hermite();

%in case we supplied bad covariates
cpp(isnan(cpp)) = 0;
cpp(isinf(cpp)) = 0; 

%index random effect categories (we expect these to be supplied in cpp consistent with X.c)
if isfield(P, 'random_effect_cat_covar_index'),
  rpp = cpp(:, P.random_effect_cat_covar_index);
  if size(cpp, 2) > 1,
    cpp(:, P.random_effect_cat_covar_index) = [];
  else
    cpp(:, P.random_effect_cat_covar_index) = 0;
  end
else
  rpp = ones(size(xpp));
end

%if we did not specify exposures, zero out ypp
if ~exist('ypp', 'var'), ypp = 0; end

tmp = cell(size(xpp));
for x = [1:length(xpp); xpp'],
  i = x(1); j = x(2);
  tmp{i} = [i*ones(j + 1, 1) (1:(j + 1))'];
end
tmp = cat(1, tmp{:});

idx = sparse(tmp(:, 1), tmp(:, 2), 1);
[pps_i pps_ct] = find(idx);
pps_i = as_column(pps_i); pps_ct = as_column(pps_ct);
pps_ct = pps_ct - 1;

[pps_i, si] = sort(pps_i);
pps_ct = pps_ct(si);
pps_midx = [diff(pps_i) == 1; true];

post_pred_surv_buf = NaN(nnz(idx), P.pp_buffer_size);
pps_max_buf = NaN(length(xpp), P.pp_buffer_size);
pps_min_buf = NaN(length(xpp), P.pp_buffer_size);
post_pred_buf = NaN(length(xpp), P.pp_buffer_size); 

%
%loop over posterior draws, averaging PP
count_PP = 0;
count_PPavg = 0;
for i = P.burnin:P.niter,
  if ~mod(i, P.thin) || i == P.niter,
    eBC = cpp*Beta(:, i) + ypp;

    buf_pos = mod(count_PP, P.pp_buffer_size) + 1;

    %compute all probabilities for CDF
    post_pred_surv_buf(:, buf_pos) = ppdraw(pps_ct, eBC(pps_i), mu(rpp(pps_i), i), 1./(1./tau(rpp(pps_i), i) + 1./tau(end, i)), X);

    %flush buffer
    if ~mod(count_PP + 1, P.pp_buffer_size) || i == P.niter,
      %cap to machine epsilon
      epsidx = post_pred_surv_buf > log(eps);

      %calculate p-values (survival functions)
      for j = 1:buf_pos,
	idx = ~pps_midx & epsidx(:, j);
	pps_max_buf(:, j) = accumarray(pps_i(idx), post_pred_surv_buf(idx, j), [length(xpp) 1], ...
	                               @(x) 1 - exp(max(x) + log(sum(exp(x - max(x))))), 1);
	pps_min_buf(:, j) = accumarray(pps_i(epsidx(:, j)), post_pred_surv_buf(epsidx(:, j), j), [length(xpp) 1], ...
	                               @(x) 1 - exp(max(x) + log(sum(exp(x - max(x))))));
	post_pred_buf(:, j) = exp(post_pred_surv_buf(pps_midx, j));
      end

      %apply machine eps. cap
      pps_max_buf(pps_max_buf <= 0) = eps;
      pps_min_buf(pps_min_buf <= 0) = eps;

      %iterative average pvalues/probabilities
      pps = cat(3, [quantile(pps_max_buf(:, 1:buf_pos), [0.025 0.975], 2) mean(pps_max_buf(:, 1:buf_pos), 2)], ...
                   [quantile(pps_min_buf(:, 1:buf_pos), [0.025 0.975], 2) mean(pps_min_buf(:, 1:buf_pos), 2)]);
      pp = [quantile(post_pred_buf(:, 1:buf_pos), [0.025 0.975], 2) mean(post_pred_buf(:, 1:buf_pos), 2)];

      if count_PPavg == 0,
	w = 1;

	pvalues = pps;
	post_pred = pp;
      else
	%if we are flushing the buffer before it's full (i.e. last flush!), we need to weight its
	%contribution.
	if buf_pos == P.pp_buffer_size, w = 1;
	else w = buf_pos/P.pp_buffer_size; end

	pvalues = pvalues + w/(count_PPavg + w)*(pps - pvalues);
	post_pred = post_pred + w/(count_PPavg + w)*(pp - post_pred);
      end

      count_PPavg = count_PPavg + w;
    end

    count_PP = count_PP + 1;

    fprintf('[%d,%d/%d] ', count_PP, i, P.niter);
  end
end

end

function pp = ppdraw(x, eBC, mu, tau, X)
  pp = -gammaln(x + 1) - log(pi)/2 + x.*(eBC + mu) + pois_LN_int(x, eBC, mu, tau, X);
end

function g = pois_LN_int(x, eBC, mu, tau, X)
  I = (sqrt(2./tau).*x)*X.h_r' - exp(bsxfun(@plus, eBC + mu, sqrt(2./tau)*X.h_r'));
  m = max(I, [], 2);
  I = I - repmat(m, 1, 100);

  g = m + log(exp(I)*X.q_w);
end
