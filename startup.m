addpath('./funcs')

[~, p] = unix('find src ref -type d');
p = strsplit(p, '\n');

for i = p, addpath(decell(i)); end
