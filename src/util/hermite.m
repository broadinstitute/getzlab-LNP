function [r, w] = hermite(file)

%loads precalcualted roots/quadrature weights for Hermite polynomials.
%can accept two column TSV or binary storing doubles

if ~exist('file'),
  file = 'hermite_roots.txt';
end

%TODO: less hacky version to determine if file is ASCII or binary.
%perhaps checking the first n bytes of file and seeing if any are > 127
if all(file((end - 2):end) == 'txt'),
  isbin = 0;
elseif all(file((end - 2):end) == 'bin'),
  isbin = 1;
else
  error('Unknown format for lookup table!');
end

if isdeployed,
  file = [ctfroot '/ref/' file];
end

if ~isbin,
  LuT = load_matrix(file);
else
  error('Binary LuT not yet implemented.');
end

r = LuT(:, 1); 
w = LuT(:, 2);
