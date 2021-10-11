%%
if(false)
mycase = 'norne'
switch mycase
    case 'norne'
        matrixfile = '/data/hnil/BITBUCKET/mrst-bitbucket/projects/project-multiscale-preformance/opm_scripts/norne_pull/reports/prob_4_4.8384e+06_matrix_istl.txt'
        rhsfile = '/data/hnil/BITBUCKET/mrst-bitbucket/projects/project-multiscale-preformance/opm_scripts/norne_pull/reports/prob_4_4.8384e+06_rhs_istl.txt'
    case 'simple'
        matrixfile = 'matrix_istl.txt'
        rhsfile = 'rhs_istl.txt'
end

mat_org = readMatrixMarket(matrixfile,true);
rhs_org = readMatrixMarket(rhsfile,false);
end
%%
mat=mat_org;
rhs=rhs_org;
d = rhs*0;
n = numel(d);
d(:)=1;
d(2:3:end)=1000*barsa;
dd = rhs*0;
dd(:)=1;
dd(3:3:end)=1/600;
Sp = sparse(1:n,1:n,d,n,n);
Sg = sparse(1:n,1:n,dd,n,n);
mat=Sg*mat_org*Sp;
rhs=rhs;
%
tol=1e-6
maxiter=40
simple_prec= struct('preconditioner','Jac','w',1,'n',1);
%amg = struct('maxlevel',5,'coarsenTarget',4000,'smoother','ILU0','alpha',0.67,'beta',1e-4);
amg = struct('maxlevel',5,'coarsenTarget',1000,'smoother','ILU0', ...
             'alpha',0.2,'beta',1e-4,'verbosity',0,'n',1,'w',1);
%Jac failed on norne
%amg_solver = struct('preconditioner','cpr','w',1.0,'n',1,'amg',amg);
%coarsesolver_simple = struct('tol',1e-2,'maxiter',1,'preconditioner','ILU0','w',1.0,'n',1)
coarsesolver_amg = struct('tol',1e-1,'maxiter',20,'preconditioner','amg','amg',amg,'verbosity',0,'solver','bicgstab');
cpr = struct('finesmoother',simple_prec,'coarsesolver',coarsesolver_amg,'verbosity',11);
opt=struct('preconditioner','cpr','w',1.0,'n',1,'amg',[],'cpr',cpr,'verbosity',10,'solver','gmres','restart',20);
tic
[x_o,extra] = duneistl(mat,rhs,'blocksize',3,'tol',tol,'maxiter',maxiter, ...
                     'istloptions',opt);
toc
%str2num(extra.time.solvetotal)
%str2num(extra.time.solvetotal)
extra.time

%%{
tic
optt=struct('preconditioner','cprt','w',1.0,'n',1,'amg',[],'cpr',cpr,'verbosity',10,'solver','gmres','restart',20);
tmat=mat';
[x_t,extra] = duneistl(tmat,rhs,'blocksize',3,'tol',tol,'maxiter',maxiter, ...
                     'istloptions',optt);
toc
%str2num(extra.time.solvetotal)
%str2num(extra.time.solvetotal)
extra.time
%}
return
norm(mat*x_o-rhs,inf)/norm(rhs,inf)
%norm(tmat*x_t-rhs)/norm(rhs)
w_cpr = readMatrixMarket("weight_cpr.txt",false);
w_cprt = readMatrixMarket("weight_cprt.txt",false);
a=reshape(w_cpr,3,[])';
ind=find(any(a<0,2));
cl = mcolon(3*(ind-1)+1,3*(ind-1)+1+2);
lcl = true(size(mat,1),1);
lcl(cl)=false;
nmat=mat(lcl,lcl);
nrhs=rhs(lcl);
[x_o,extra] = duneistl(nmat,nrhs,'blocksize',3,'tol',tol,'maxiter',maxiter, ...
                       ...
'istloptions',opt);
%[x_t,extra] = duneistl(nmat',nrhs,'blocksize',3,'tol',tol,'maxiter',maxiter, ...
%                       ...
%                     'istloptions',optt);

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
