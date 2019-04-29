function x = duneistl(mat, rhs, varargin)
opt = struct('blocksize',1,...
             'tol', 1e-3,...
             'maxiter',200);
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
x = duneistl_matlab(i,j,val, rhs, opt.blocksize, opt.tol, opt.maxiter); 
end