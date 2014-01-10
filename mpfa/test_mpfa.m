function test_mpfa()
%Example demonstrating basic use of the MPFA-O pressure solver.

require mimetic

gravity off

g = cartGrid([30, 30]);
g.nodes.coords = twister(g.nodes.coords);
g = computeGeometry(g);
%--------------------------------------------------------------------------

clear rock
for i = 1:size(g.nodes.coords, 2),
   rock.perm(:,i) = convertFrom(logNormLayers([g.cartDims, 1], 1), ...
                                milli*darcy);
end
rock = struct('perm', repmat(0.1*darcy, [g.cells.num, 1]));

%--------------------------------------------------------------------------

%{
src = addSource([], [1 ; g.cells.num], [1, -1]*meter^3/second);
%}
%%{
src = [];
%}

bc  = pside([], g, 'left',  2);
bc  = pside(bc, g, 'right', 1);

%--------------------------------------------------------------------------

fprintf('Mimetic Method\t... ')
tic
S = computeMimeticIP(g, rock);
fluid = initSingleFluid('mu' ,    1*centi*poise     , ...
                        'rho', 1014*kilogram/meter^3);
xr1 = solveIncompFlow(initResSol(g, 0, 0), g, S, fluid, ...
                      'src', src, 'bc', bc);
toc

%--------------------------------------------------------------------------

fprintf('MPFA O-Method\t... ')
tic
T1  = computeMultiPointTrans(g, rock);
xr2 = incompMPFA(initResSol(g, 0, 0), g, T1, fluid, ...
                 'src', src, 'bc', bc,'MatrixOutput',true);
toc

%--------------------------------------------------------------------------

fprintf('TPFA Method\t... ')
tic
T2  = computeTrans(g, rock);
xr3 = incompTPFA(initResSol(g, 0, 0), g, T2, fluid, ...
                 'src', src, 'bc', bc,'MatrixOutput',true);
toc

%--------------------------------------------------------------------------

cellno     = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
flux_int   = @(x) accumarray(cellno, abs(x.flux(g.cells.faces(:,1))));

plot_var   = @(v) plotCellData(g, v);
%{
X          = reshape(g.cells.centroids(:,1), g.cartDims(1:2));
Y          = reshape(g.cells.centroids(:,2), g.cartDims(1:2));
plot_press = @(x) contourf(X, Y, ...
                           reshape(x.cellPressure, g.cartDims(1:2)), ...
                           linspace(min(bc.value), max(bc.value), 11));
%}
%%{
plot_press = @(x) plot_var(x.pressure(1:g.cells.num));
%}
plot_flux  = @(x) plot_var(convertTo(flux_int(x), meter^3/day));

clf
subplot(2,3,1),
   plot_flux(xr1); colorbar, cax = caxis;
   axis equal tight, title('Mimetic')

subplot(2,3,2),
   plot_flux(xr2); colorbar, caxis(cax);
   axis equal tight, title('MPFA-O')

subplot(2,3,3),
   plot_flux(xr3); colorbar, caxis(cax);
   axis equal tight, title('TPFA')

subplot(2,3,4),
   plot_press(xr1); colorbar, cax = caxis;
   axis equal tight, title('Mimetic')

subplot(2,3,5),
   plot_press(xr2); colorbar, caxis(cax);
   axis equal tight, title('MPFA-O')

subplot(2,3,6),
   plot_press(xr3); colorbar, caxis(cax);
   axis equal tight, title('TPFA')

%--------------------------------------------------------------------------

p.pressure = 2 - g.cells.centroids(:,1)/g.cartDims(1);

err        = @(q1, q2) norm(q1 - q2, inf);
err_press  = @(x1, x2) err(x1.pressure(1:g.cells.num), ...
                           x2.pressure(1:g.cells.num));
err_flux   = @(x1, x2) err(flux_int(x1), flux_int(x2));

fprintf(['\nFlux Difference:\n', ...
         '\to Mimetic/MPFA-O : %.15e\n',    ...
         '\to Mimetic/TPFA   : %.15e\n',    ...
         '\to MPFA-O /TPFA   : %.15e\n\n'], ...
        err_flux(xr1, xr2), err_flux(xr1, xr3), err_flux(xr2, xr3));

fprintf(['Cell Pressure Difference:\n', ...
         '\to Mimetic/MPFA-O : %.15e\n',    ...
         '\to Mimetic/TPFA   : %.15e\n',    ...
         '\to MPFA-O /TPFA   : %.15e\n\n'], ...
        err_press(xr1, xr2), err_press(xr1, xr3), err_press(xr2, xr3));

%%{
fprintf(['Cell Pressure Error:\n', ...
         '\to Mimetic        : %.15e\n',  ...
         '\to MPFA-O         : %.15e\n',  ...
         '\to TPFA           : %.15e\n'], ...
        err_press(xr1, p), err_press(xr2, p), err_press(xr3, p));
%}
