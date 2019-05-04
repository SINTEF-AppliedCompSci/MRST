function [x,extra] = duneistl(mat, rhs, varargin)
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
    simple_prec= struct('preconditioner','ILU0','w',1.0,'n',1);
    amg = struct('maxlevel',4,'coarsenTarget',1000,'smoother','ILU0','alpha',0.67,'beta',1e-4);
    amg_solver = struct('preconditioner','cpr','w',1.0,'n',1,'amg',amg);
    coarsesolver_simple = struct('tol',1e-2,'maxiter',100,'preconditioner','ILU0','w',1.0,'n',1)
    coarsesolver_amg = struct('tol',1e-2,'maxiter',100,'preconditioner','amg','amg',amg);
    cpr = struct('finesmoother',simple_prec,'coarsesolver',coarsesolver_amg);
    sopt=struct('preconditioner','cpr','w',1.0,'n',1,'amg',amg,'cpr',cpr); 
    options = jsonencode(sopt);
else
    options = jsonencode(opt.istloptions);
end

[x,ext] = duneistl_matlab(i,j,val, rhs, opt.blocksize, opt.tol, ...
                          opt.maxiter, options);
extra = jsondecode(ext);
end