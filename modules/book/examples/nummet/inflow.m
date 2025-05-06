function u=inflow(u,n)
if nargin==1
   u(1) = 1;
else
   u(1:n) = 1;
end