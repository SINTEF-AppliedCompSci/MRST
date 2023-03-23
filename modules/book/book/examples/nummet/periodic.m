function u = periodic(u,n)
if nargin==1
    u(1) = u(end-1);
    u(end) = u(2);
else
    u(1:n) = u(end-2*n+1:end-n);
    u(end-n+1:end) = u(n+1:2*n);
end