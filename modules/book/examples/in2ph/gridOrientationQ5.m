%% Grid-orientation effects for the five-spot
% With a two-point type discretization (and with many other discretizations
% as well), evolving displacement profiles will preferrentially move along
% the axial directions, i.e., in the direction of the normals to the cell
% faces. To illustrate this, we contrast approximate solutions for the
% repeated five-spot well pattern computed using the original quarter
% five-spot setup and with a rotated setup with injectors in the SW and NE
% corners and producers in the SE and NW corners. The two setups will have
% different preferrential flow directions and hence generally give
% different solutions.
mrstModule add incomp

%% Difference in grid-orientation effects with mobility ratio
% This effect is typically associated with instable displacements. To show
% this we consider three different fluid models with unfavorable mobility
% ratio, equal viscosities, and favorable mobility ratio.
mu       = [1 10; 1 1; 10 1];
pvi      = [0.3 0.6 0.7];
cartDim  = [32 32 1];
nstep    = 16;
set(gcf,'Position',[250 450 1100 360]);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).*centi*poise, ...
        'rho', [1000,  850].* kilogram/meter^3, 'n', [2, 2]);
    subplot(1,3,n)
    runQ5DiagParal(cartDim, nstep, fluid, pvi(n));
    set(gca,'XTick',[],'YTick',[]);
    title(sprintf('Mobility ratio %d:%d. Time: %.1f PVI',...
        mu(n,1),mu(n,2),pvi(n)),'FontSize',12); drawnow
end

%% Convergence with respect to time step
m = 1;
fluid = initSimpleFluid('mu', mu(m,:).*centi*poise, ...
    'rho', [1000,  850].* kilogram/meter^3, 'n', [2, 2]);
nsteps = [1 8 64];
for n=1:3
     subplot(1,3,n)
     runQ5DiagParal(cartDim, nsteps(n), fluid, pvi(m));
     set(gca,'XTick',[],'YTick',[]);
     title(sprintf('Grid: %d x %d. Steps: %d', ...
         cartDim(1:2), nsteps(n)),'FontSize',12); drawnow
end

%% Convergence with respect to spatial resolution
m = 1;
fluid = initSimpleFluid('mu', mu(m,:).*centi*poise, ...
    'rho', [1000,  850].* kilogram/meter^3, 'n', [2, 2]);
cartDims = [16 16 1; 32 32 1; 64 64 1];
for n=1:3
     subplot(1,3,n)
     runQ5DiagParal(cartDims(n,:), nstep, fluid, pvi(m));
     set(gca,'XTick',[],'YTick',[]);
     title(sprintf('Grid: %d x %d. Steps: %d', ...
         cartDims(n,1:2), nstep),'FontSize',12); drawnow
end
