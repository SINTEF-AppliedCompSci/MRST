function [x,extra] = duneistl(mat, rhs, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt = struct('blocksize',1,...
             'tol', 1e-3,...
             'maxiter',200,...
              'istloptions',[]);
[opt, cl_opts] = merge_options(opt, varargin{:});

[i,j,val] = find(mat);
i=i-1;
j=j-1;

ms = size(mat)
vs =size(rhs);
assert(ms(1)==ms(2));
assert(ms(1)==vs(1));
assert(vs(2)==1);
assert(mod(vs(1),opt.blocksize) == 0);
if(isempty(opt.istloptions))
    simple_prec= struct('preconditioner','ILU0','w',1,'n',1);
    amg = struct('maxlevel',5,'coarsenTarget',1000,...
                 'smoother','ILU0', ...
                 'alpha',0.2,'beta',1e-4,'verbosity',0,'n',1,'w',1);
    %Jac failed on norne
    coarsesolver_amg = struct('tol',1e-1,'maxiter',20, ...
                              'preconditioner','amg',...
                              'amg',amg,...
                             'verbosity',0,...
                              'solver','bicgstab','pressure_var_index',1);
    cpr = struct('finesmoother',simple_prec,'coarsesolver', ...
                 coarsesolver_amg,'verbosity',11, 'pressure_var_index',1);
    sopt=struct('preconditioner','cpr','w',1.0,'n',1,'amg',[],'cpr',cpr,'verbosity',10,'solver','gmres','restart',20);
    options = jsonencode(sopt);
else
    options = jsonencode(opt.istloptions);
end
bb=opt.blocksize;
rows=[floor(i/bb),floor(j/bb),mod(i,bb),mod(j,bb)];
[~,ind] = sortrows(rows);

%[i,ind] = sort(i);
%%{
i=i(ind);
j=j(ind);
val=val(ind);
%}
[x,ext] = duneistl_matlab(i,j,val, rhs, opt.blocksize, opt.tol, ...
                          opt.maxiter, options);
extra = jsondecode(ext);
end
