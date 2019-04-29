addpath('/data/hnil/BITBUCKET/mrst-bitbucket/projects/opm-challenge/scripts/')
%%
delete('duneistl_matlab_file.mexa64');

[rhs, bzv, bmat] = readMatrixMarket('rhs_istl.txt',false);
m = numel(rhs);
bz=bzv(1);

a=duneistl_matlab_file('matrix_istl.txt','rhs_istl.txt', m, bz)
%%
[mat, bzm, bmat] = readMatrixMarket('matrix_istl.txt',true);
writeMatrixMarket(mat, bzm, 'matrix.txt',true)

writeMatrixMarket(rhs, bzv, 'rhs.txt',false)

[mat2, bz2] = readMatrixMarket('matrix.txt',true);
%%
% ISTL_STRUCT blocked 3 3

mat = readMatrixMarket('matrix_istl.txt',true)
rhs = readMatrixMarket('rhs_istl.txt',false)

am = mat\rhs-a
a2=duneistl_matlab_file('matrix.txt','rhs.txt', m, bz)

writeMatrixMarket(mat, [1,1], 'matrix_1b.txt',true);
writeMatrixMarket(rhs, [1 1], 'rhs_1b.txt',false);
a_1b=duneistl_matlab_file('matrix_1b.txt','rhs_1b.txt', m, 1)
%%
[i,j,val] = find(mat);
x=sparse(i,j,val, size(rhs,1), size(rhs,1))\rhs;    
i=i-1;
j=j-1;
%delete('duneistl_matlab.mexa64');
a_g = duneistl_matlab(i,j,val, rhs, 1, 1e-4, 200); 
