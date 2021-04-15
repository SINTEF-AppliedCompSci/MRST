addpath('/data/hnil/BITBUCKET/mrst-bitbucket/projects/opm-challenge/scripts/')
%%
%delete('duneistl_matlab_file.mexa64');

[rhs, bzv, bmat] = readMatrixMarket('rhs_istl.txt',false);
m = numel(rhs);
bz=bzv(1);
opt=struct('preconditioner','ILU0','w',1.0,'n',1,'solver','bicgstab','verbosity',10);
opts = jsonencode(opt);
a=duneistl_matlab_file('matrix_istl.txt','rhs_istl.txt', m, bz,opts)
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
a2=duneistl_matlab_file('matrix.txt','rhs.txt', m, bz, opts)

writeMatrixMarket(mat, [1,1], 'matrix_1b.txt',true);
writeMatrixMarket(rhs, [1 1], 'rhs_1b.txt',false);
a_1b=duneistl_matlab_file('matrix_1b.txt','rhs_1b.txt', m, 1,opts)

% [i,j,val] = find(mat);
% x=sparse(i,j,val, size(rhs,1), size(rhs,1))\rhs;    
% i=i-1;
% j=j-1;
%delete('duneistl_matlab.mexa64');
%a_g = duneistl_matlab(i,j,val, rhs, 1, 1e-4, 200); 
%%
simple_prec= struct('preconditioner','ILU0','w',1.0,'n',1,'verbosity',10);
amg = struct('maxlevel',4,'coarsenTarget',1000,'smoother','ILU0','alpha',0.67,'beta',1e-4,'verbosity',10,'n',1,'w',1);
amg_solver = struct('preconditioner','cpr','w',1.0,'n',1,'amg',amg,'verbosity',10);
coarsesolver_simple = struct('tol',1e-2,'maxiter',100,'preconditioner','ILU0','w',1.0,'n',1,'verbosity',10)
coarsesolver_amg = struct('tol',1e-2,'maxiter',100,'preconditioner','amg','amg',amg,'verbosity',10,'n',1,'w',1,'solver','bicgstab');
cpr = struct('finesmoother',simple_prec,'coarsesolver',coarsesolver_amg,'verbosity',10);
opt=struct('preconditioner','cpr','w',1.0,'n',1,'amg',amg,'cpr',cpr,'verbosity',10,'solver','bicgstab','restart',10);   
x = duneistl(mat,rhs,'blocksize',3,'tol',1e-4,'maxiter',100,'istloptions',opt)

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
