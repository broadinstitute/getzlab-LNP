function p = normgampdf(mu, tau, A, B, M, T, logmode)

%argument checks
if ~isequal(size(mu), size(tau)),
  error('Size of mu/tau inputs must be the same!');
end
if ~isequal(size(A), size(B), size(M), size(T)),
  error('Size of A/B/M/T parameters must be the same!');
end

if ~exist('logmode'), logmode = false; end
logmode = logical(logmode);

if ~logmode,
  p = tau.^(A - 1).*exp(-tau./B).*exp(-T.*tau/2.*(mu - M).^2).*sqrt(tau);
  z = B.^A.*sqrt(2*pi).*gamma(A)./sqrt(T);

  p = p/z;
else
  p = (A - 1).*log(tau) - tau./B - T.*tau/2.*(mu - M).^2 + log(tau)/2;
  z = A.*log(B) + log(sqrt(2*pi)) + gammaln(A) - log(T)/2;

  p = p - z;
end
