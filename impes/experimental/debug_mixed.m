%%
OP = tpfaConnectionOperators(g, W, 2);
x1 = x;
X  = [g.cells.centroids(:,1); g.faces.centroids(:,1)];
it = 0;

f1 = [x1.flux ; vertcat(x1.wellSol.flux)];
p1 = [x1.pressure ; vertcat(x1.wellSol.pressure)];

%%
f0 = f1;
p0 = p1;
[A, b, dofPos, fmob, rho] = impesTPFAMixedSystem(x1, g, T, fluid, DT, ...
                                                 PV, 'wells', W);

soln = A \ b;

x1 = impesMixedUpdatePressure(x1, soln, dofPos, g, T, ...
                              fmob, rho, OP, W, []);

f1 = [x1.flux ; vertcat(x1.wellSol.flux)];
p1 = [x1.pressure ; vertcat(x1.wellSol.pressure)];

dp = p1 - p0;
df = f1 - f0;

it = it + 1;

fmt = ['  %02d  ', repmat('%22.15e ', [1, 4]), '\n'];
fprintf(fmt, it, norm(dp,inf), norm(dp,inf)./norm(p0,inf), ...
                 norm(df,inf), norm(df,inf)./norm(f0,inf));
%%
Y = convertTo([x1.pressure; x1.facePressure], barsa);
V = sortrows([X, Y]);

plot(V(:,1), V(:,2), '-o');

%%
figure
plot(g.faces.centroids(:,1), x1.flux, '-*');
