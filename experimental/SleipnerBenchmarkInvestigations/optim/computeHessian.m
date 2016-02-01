function H = computeHessian(f, u, dims, abspert)
% compute sub-Hessian of f for dims
if nargin < 4
    relpert = 1e-5;
    abspert = abs(u(dims))*relpert;    
end

n = numel(dims);
H = zeros(n,n);

[v, g] = f(u);
g = g(dims);
for k = 1:numel(dims);
    d = zeros(size(u));
    d(dims(k)) = abspert(k);
    [vk, gk] = f(u+d);
    H(:,k) = (gk(dims)-g)/abspert(k);
end
end
    

   

