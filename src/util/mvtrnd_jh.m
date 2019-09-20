function x = mvtrnd_jh(mu, sigma, nu, cases, skip)

if ~exist('cases', 'var'), cases = 1; end
if ~exist('skip', 'var'), skip = false; end

if ~skip,
  %sanity checks
  if ndims(sigma) ~= 2,
    error('Covariance matrix must be a 2D!')
  end

  sz = size(sigma); 
  if sz(1) ~= sz(2),
    error('Covariance matrix must be square!')
  end

  if ~isequal([sz(1) 1], size(mu)),
    error('Mean must be a column vector with matching dimensions!')
  end

  if ~issymmetric(sigma) || ~all(eig(sigma) >= 0),
    error('Covariance matrix must be symmetric and positive-semidefinite!');
  end
end

m = size(sigma, 1);

[V l] = eig(sigma);

x = mu + V*sqrt(l)*trnd(nu, m, cases);
