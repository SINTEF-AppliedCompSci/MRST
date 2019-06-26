function [x,extra] = amgcl(mat, rhs, varargin)
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
    options = struct('solver',solver,'precond',precond,'solver_type','regular','write_params',true,'block_size',opt.blocksize,'verbose',true);    
else
    options = jsonencode(opt.istloptions);
end
options = jsonencode(options);
[x, err, nIter] = amgcl_matlab_simple(mat', rhs, opt.tol, opt.maxiter, options);
extra =struct('err',err,'nIter',nIter);                        
