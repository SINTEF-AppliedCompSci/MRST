clear all

test = [struct('index', [1 2 3]), struct('index', [4 5])];

[a, b] = initVariablesADI(ones(10, 1), zeros(2,1));

%%
clc
% Should throw error
disp('Multiindex')
a(test.index)
%%
clc
% Should be ok
disp('Subindex of index')
a(1).jac

%%
clc
% Should be ok
disp('Single index')
a(1)
