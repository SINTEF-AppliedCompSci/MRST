% Loading specifications of system
system_specification;


T = 1*day;
r = 1:0.1:10000;

% Constructing analytical functions
[pfun, hfun, zones_fun] = ...
    similarity_sol(H, phi, k, Qvol, Sres, krc, muc, mub, rhoc, rhob, c_r, c_b);

% Testing
p = pfun(r,T,10000,T, 0) + Bpress;
%p = hfun(r,T);
plot(log10(r), p);
hold on;
zones = zones_fun(T);
for i = 1:3
    plot(log10([zones(i), zones(i)]), [min(p), max(p)], 'r');
end
zones_fun(T)
hold off;
%plot(r, p);