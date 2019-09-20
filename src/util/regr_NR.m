function [mu, tau, Beta, H] = regr_NR(X, Y, C, HP, P)

if nargin < 3,
  error('Too few arguments!')
elseif nargin <= 4, %no exposures
  C = Y;
  HP = C;

  if exist('HP', 'var'), P = HP; end
elseif nargin >= 5,
  %nothing to do!
end

%performs Newton-Raphson or gradient descent iterations to find MAP for regression model

C(isnan(C)) = 0;
C(isinf(C)) = 0;

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

%ensure hyperparameters are present
%                  normal-gamma                 multivariate normal
demand_fields(HP, {'taua' 'taub' 'mumu' 'mutau' 'invcovmu' 'mubeta'})
HP = rename_fields(HP, {'taua' 'taub' 'mumu' 'mutau' 'invcovmu' 'mubeta'}, ...
                       {'A' 'B' 'M' 'T' 'invcovmat' 'muvec'});

P = impose_default_value(P, 'method', 'NR');
P = impose_default_value(P, 'min_gradient_norm', 1e-8);
P = impose_default_value(P, 'learnrate', 2e-5);
P = impose_default_value(P, 'NRlearnrate', 1);
if strcmp(P.method, 'GD'),
  P = impose_default_value(P, 'niter', 1000);
elseif strcmp(P.method, 'NR'), 
  P = impose_default_value(P, 'niter', 20);
else
  error(['Invalid method: ' P.method]);
end
P = impose_default_value(P, 'm0', -2);
P = impose_default_value(P, 't0', 1);
P = impose_default_value(P, 'beta0', zeros(size(C, 2), 1));

[h_r, q_w] = hermite();
q_wD = sparse(1:length(q_w), 1:length(q_w), q_w);

mu = NaN(P.niter, 1); mu(1) = P.m0;
tau = NaN(P.niter, 1); tau(1) = log(P.t0);
Beta = NaN(P.niter, size(C, 2)); Beta(1, :) = P.beta0;
i = 2;
while true,
  %
  %compute common terms
  tau_e = exp(tau(i - 1));

  eBC = exp(C*Beta(i - 1, :)' + Y);
  logI = qvec(X, eBC, mu(i - 1), tau_e, h_r);

  md = max(logI, [], 2);
  logI = bsxfun(@minus, logI, md);

  %note that the term pulled in logsumexp will always cancel, as sums are always part of quotients.
  I_mat = exp(logI)*q_wD;
  I_vec = sum(I_mat, 2);

  %
  %compute gradient 
  g = zeros(size(C, 2) + 2, 1);
  g(1) = gmu(I_mat, I_vec, mu(i - 1), tau_e, h_r, HP);
  g(2) = gtau(I_mat, I_vec, mu(i - 1), tau_e, length(X), h_r, HP);
  g(3:end) = gbeta(I_mat, I_vec, mu(i - 1), tau_e, Beta(i - 1, :)', eBC, C, X, h_r, HP);

  if strcmp(P.method, 'NR'), 
    %
    %compute Hessian elements
    H = zeros(size(C, 2) + 2);

    %off-diagonal elements
    H(2, 1) = hmutau(I_mat, I_vec, mu(i - 1), tau_e, h_r, HP)*tau_e;

    for j = 3:(size(C, 2) + 2),
      H(j, 1) = hmubeta(I_mat, I_vec, mu(i - 1), tau_e, eBC, C(:, j - 2), X, h_r);
    end

    for j = 3:(size(C, 2) + 2),
      H(j, 2) = htaubeta(I_mat, I_vec, mu(i - 1), tau_e, eBC, C(:, j - 2), X, h_r)*tau_e;
    end

    for j = 3:(size(C, 2) + 2), for k = (j + 1):(size(C, 2) + 2),
      H(k, j) = hbetabeta(I_mat, I_vec, tau_e, eBC, C(:, j - 2), C(:, k - 2), X, h_r, HP.invcovmat(k - 2, j - 2));
    end, end

    H = H + H';

    %on-diagonal elements
    H(1, 1) = hmumu(I_mat, I_vec, mu(i - 1), tau_e, h_r, HP);
    H(2, 2) = htautau(I_mat, I_vec, mu(i - 1), tau_e, length(X), h_r, HP)*tau_e^2 + g(2)*tau_e;

    for j = 3:(size(C, 2) + 2),
      H(j, j) = hbetabeta(I_mat, I_vec, tau_e, eBC, C(:, j - 2), C(:, j - 2), X, h_r, HP.invcovmat(j - 2, j - 2));
    end

    %
    %change-of-variable chain rule factor for gradient
    g(2) = g(2)*tau_e; g(2) = g(2) + 1;

    %
    %update parameters
    step = inv(H)*g;

    mu(i) = mu(i - 1) - P.NRlearnrate*step(1);
    tau(i) = tau(i - 1) - P.NRlearnrate*step(2);
    Beta(i, :) = Beta(i - 1, :) - P.NRlearnrate*step(3:end)';
  elseif strcmp(P.method, 'GD'),
    mu(i) = mu(i - 1) + P.learnrate*g(1);
    tau(i) = tau(i - 1) + P.learnrate*g(2);
    Beta(i, :) = Beta(i - 1, :) + P.learnrate*g(3:end)';
  end

  %
  %summarize iteration
  fprintf('post = %0.4f ||g|| = %0.2e\n', ...
          -length(X)/2*log(pi) + ...
          log(eBC)'*X - length(X)*tau_e*mu(i - 1)^2/2 - sum(gammaln(X + 1)) + ...
          sum(md + log(I_vec)) + ...
          normgampdf(mu(i - 1), tau_e, HP.A, HP.B, HP.M, HP.T, true) + ...
          log(mvnpdf(Beta(i - 1, :)', HP.muvec, inv(HP.invcovmat))), norm(g));

  fprintf('%0.4f ', [exp(tau(i)) mu(i)]);
  fprintf('\n');

  fprintf('%0.2f ', Beta(i, :));
  fprintf('\n\n');

  %
  %stop if we've reached maximum number of iterations or a maximum
  if i >= P.niter || (norm(g) <= P.min_gradient_norm && all(eig(H) < 0)),
    break;
  end

  i = i + 1;
end

mu = mu(1:i);
tau = exp(tau(1:i));
Beta = Beta(1:i, :);

end

%
% GRADIENT FUNCTIONS {{{
%

%
%mu
function g = gmu(I_mat, I_vec, mu, tau, h_r, H)
  g = sum((I_mat*hB(mu, tau, h_r))./I_vec) ...
      - H.T*tau*(mu - H.M);
end

% 
% tau
function g = gtau(I_mat, I_vec, mu, tau, len, h_r, H)
  g = sum((I_mat*-hA(mu, tau, h_r))./I_vec) + len/(2*tau) ...
      + (H.A - 1)/tau - 1/H.B - H.T/2*(mu - H.M)^2 + 1/(2*tau);
end

%
% beta
function g = gbeta(I_mat, I_vec, mu, tau, Beta, eBC, C, X, h_r, H)
  g = -H.invcovmat*(Beta - H.muvec);
  for j = 1:size(C, 2),
    g(j) = g(j) + sum(sum(I_mat.*hC(tau, C(:, j), X, eBC, h_r), 2)./I_vec);
  end
end

% }}}

%
% HESSIAN FUNCTIONS {{{
%

%
%mu, mu 
function h = hmumu(I_mat, I_vec, mu, tau, h_r, H)
  h = sum(I_mat*(hB(mu, tau, h_r).^2 - tau)./I_vec - ((I_mat*hB(mu, tau, h_r))./I_vec).^2) ...
      - H.T*tau;
end

%
%tau, tau
function h = htautau(I_mat, I_vec, mu, tau, len, h_r, H)
  h = sum(I_mat*hA(mu, tau, h_r).^2./I_vec - ((I_mat*hA(mu, tau, h_r))./I_vec).^2) - len/(2*tau^2) ...
      - (H.A - 1)/tau^2 - 1/(2*tau^2);
end

%
%beta(j), beta(k)
function h = hbetabeta(I_mat, I_vec, tau, eBC, C1, C2, X, h_r, HPicm)
  h = sum(sum(I_mat.*(hC(tau, C1, X, eBC, h_r).*hC(tau, C2, X, eBC, h_r) - ...
              C1.*C2.*eBC*exp(sqrt(2/tau)*h_r)'), 2)./I_vec - ...
          (sum(I_mat.*hC(tau, C1, X, eBC, h_r), 2).*sum(I_mat.*hC(tau, C2, X, eBC, h_r), 2))./I_vec.^2) ...
      - HPicm;
end

%
%mu, tau
function h = hmutau(I_mat, I_vec, mu, tau, h_r, H)
  h = sum(I_mat*((1/tau - hA(mu, tau, h_r)).*hB(mu, tau, h_r))./I_vec - ...
          ((I_mat*hA(mu, tau, h_r)).*(I_mat*-hB(mu, tau, h_r)))./I_vec.^2) ...
      - H.T*(mu - H.M);
end

%
%mu, beta
function h = hmubeta(I_mat, I_vec, mu, tau, eBC, C, X, h_r)
  h = sum(sum(I_mat.*(bsxfun(@times, hB(mu, tau, h_r)', hC(tau, C, X, eBC, h_r))), 2)./I_vec - ...
          ((I_mat*hB(mu, tau, h_r)).*(sum(I_mat.*hC(tau, C, X, eBC, h_r), 2)))./I_vec.^2);
end

%
%tau, beta
function h = htaubeta(I_mat, I_vec, mu, tau, eBC, C, X, h_r)
  h = sum(sum(I_mat.*(bsxfun(@times, -hA(mu, tau, h_r)', hC(tau, C, X, eBC, h_r))), 2)./I_vec - ...
          ((I_mat*-hA(mu, tau, h_r)).*(sum(I_mat.*hC(tau, C, X, eBC, h_r), 2)))./I_vec.^2);
end 

% }}}

%
% MISC. FUNCTIONS {{{
%

%
%log quadrature vectors
function q = qvec(X, eBC, mu, tau, h_r)
  q = (sqrt(2/tau)*X + mu*sqrt(2*tau))*h_r' - eBC*exp(sqrt(2/tau)*h_r)'; 
end

%
% "Eigenvalue" functions {{{

function a = hA(mu, tau, h_r)
  a = 1/2*(2/tau*h_r.^2 - 2*mu*h_r*sqrt(2/tau) + mu^2);
end

function b = hB(mu, tau, h_r)
  b = tau*(h_r*sqrt(2/tau) - mu);
end

function c = hC(tau, C, X, eBC, h_r)
  c = bsxfun(@minus, C.*X, C.*eBC*exp(sqrt(2/tau)*h_r)');
end

% }}}

% }}}
