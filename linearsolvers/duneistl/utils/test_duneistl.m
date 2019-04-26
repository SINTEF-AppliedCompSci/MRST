addpath('/data/hnil/BITBUCKET/mrst-bitbucket/projects/opm-challenge/scripts/')
%%
delete('duneistl_matlab.mexa64');
a=duneistl_matlab('matrix_istl.txt','rhs_istl.txt')
%%
[mat, bzm, bmat] = readMatrixMarket('matrix_istl.txt',true);
writeMatrixMarket(mat, bzm, 'matrix.txt',true)
[rhs, bzv, bmat] = readMatrixMarket('rhs_istl.txt',false);
writeMatrixMarket(rhs, bzv, 'rhs.txt',false)

[mat2, bz2] = readMatrixMarket('matrix.txt',true);
%%
% ISTL_STRUCT blocked 3 3

mat = readMatrixMarket('matrix_istl.txt',true)
rhs = readMatrixMarket('rhs_istl.txt',false)

am = mat\rhs-a
a2=duneistl_matlab('matrix.txt','rhs.txt')