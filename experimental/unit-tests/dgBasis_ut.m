tol = 1e-9;

%%
degree = 10;

nDof = polyDim(degree, 1);
[l, d_l] = deal(cell(nDof,1));

l{1} = @(x) 0*x + 1;
d_l{1} = @(x) 0*x;
if degree > 0
    l{2} = @(x) x;
    d_l{2} = @(x) 0*x + 1;
    for k = 1:nDof-2
        l{k+2} = @(x) ((2*k+1)*x.*l{k+1}(x) - k*l{k}(x))./(k+1);
        d_l{k+2} = @(x) (k+1)*(x.*l{k+2}(x) - l{k+2-1}(x))./(x.^2-1);
    end
end

%%

fprintf('Testing 1D Legenrdre polynomials ... ')

dim = 1;
basis    = dgBasis(dim, degree, 'legendre');
psi      = basis.psi;
grad_psi = basis.grad_psi;
nDof     = basis.nDof;
k        = basis.k;

xl       = 0.9999;
x = linspace(-xl, xl, 100)';
fun_diff = cellfun(@(l1,l2) norm(l1(x) - l2(x)), l, psi);
grad_fun_diff = cellfun(@(gl1,gl2) norm(gl1(x) - gl2(x)), d_l, grad_psi);

assert(all(fun_diff<tol) & all(grad_fun_diff < tol));

fprintf('1D unit test successful!\n')

%%

fprintf('Testing 2D Legenrdre polynomials ... ')

dim      = 2;
basis    = dgBasis(dim, degree, 'legendre');
psi      = basis.psi;
grad_psi = basis.grad_psi;
nDof     = basis.nDof;
k        = basis.k;

[L, dL] = deal(cell(nDof,1));
for dofNo = 1:nDof
    L{dofNo} = @(x) 0*x + 1;
    for dNo = 1:dim
        L{dofNo} = @(x) L{dofNo}(x).*l{k(dofNo, dNo) + 1}(x(:,dNo));
    end
    dL{dofNo} = @(x) [d_l{k(dofNo,1)+1}(x(:,1)).*  l{k(dofNo,2)+1}(x(:,2)), ...
                        l{k(dofNo,1)+1}(x(:,1)).*d_l{k(dofNo,2)+1}(x(:,2))];
end

x = linspace(-xl, xl, 10)';
[xx,yy] = ndgrid(x, x);
x = [xx(:), yy(:)];

fun_diff = cellfun(@(l1,l2) norm(l1(x) - l2(x)), L, psi);
grad_fun_diff = cellfun(@(gl1,gl2) norm(sqrt(sum((gl1(x) - gl2(x)).^2,2))), dL, grad_psi);
assert(all(fun_diff<tol) & all(grad_fun_diff < tol));

fprintf('2D unit test successful!\n')

%%

fprintf('Testing 3D Legenrdre polynomials ... ')

dim      = 3;
basis    = dgBasis(dim, degree, 'legendre');
psi      = basis.psi;
grad_psi = basis.grad_psi;
nDof     = basis.nDof;
k        = basis.k;

[L, dL] = deal(cell(nDof,1));
for dofNo = 1:nDof
    L{dofNo} = @(x) 0*x + 1;
    for dNo = 1:dim
        L{dofNo} = @(x) L{dofNo}(x).*l{k(dofNo, dNo) + 1}(x(:,dNo));
    end
    dL{dofNo} = @(x) [d_l{k(dofNo,1)+1}(x(:,1)).*  l{k(dofNo,2)+1}(x(:,2)).*  l{k(dofNo,3)+1}(x(:,3)), ...
                        l{k(dofNo,1)+1}(x(:,1)).*d_l{k(dofNo,2)+1}(x(:,2)).*  l{k(dofNo,3)+1}(x(:,3)), ...
                        l{k(dofNo,1)+1}(x(:,1)).*  l{k(dofNo,2)+1}(x(:,2)).*d_l{k(dofNo,3)+1}(x(:,3))];
end

x = linspace(-xl, xl, 10)';
[xx,yy, zz] = ndgrid(x, x, x);
x = [xx(:), yy(:), zz(:)];

fun_diff = cellfun(@(l1,l2) norm(l1(x) - l2(x)), L, psi);
grad_fun_diff = cellfun(@(gl1,gl2) norm(sqrt(sum((gl1(x) - gl2(x)).^2,2))), dL, grad_psi);
assert(all(fun_diff<tol) & all(grad_fun_diff < tol));

fprintf('3D unit test successful!\n')