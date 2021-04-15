function [x,extra] = amgcl(mat, rhs, varargin)
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
             'amgcloptions',[]);
[opt, cl_opts] = merge_options(opt, varargin{:});
ms = size(mat)
vs =size(rhs);
assert(ms(1)==ms(2));
assert(ms(1)==vs(1));
assert(vs(2)==1);
assert(mod(vs(1),opt.blocksize) == 0);
if(isempty(opt.amgcloptions))
    solver = struct('type','gmres','M',20);
    precond = struct('class','relaxation','type','ilu0','damping',1);
    options = struct('solver',solver,'precond',precond,'solver_type','regular',...
        'write_params',true,'block_size',opt.blocksize,'verbosity',0,'reuse_mode',1);    
else
    options = opt.amgcloptions;
end
options = jsonencode(options);
[x, err, nIter] = amgcl_matlab_simple(mat', rhs, opt.tol, opt.maxiter, options);
extra =struct('err',err,'nIter',nIter);                        
