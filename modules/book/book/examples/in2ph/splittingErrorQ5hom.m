%% Quarter five-spot: Illustration of temporal splitting errors
% We use the standard quarter five-spot test case to illustrate splitting
% errors that may arise when using a sequential solution procedure to
% simulate multiphase flow 
mrstModule add book incomp diagnostics

%% Set up model
T = 1;
cartDim = [128 128 1];
gravity reset off
fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
domain = [250 250 20];
G      = computeGeometry(cartGrid(cartDim,domain));
rock   = makeRock(G, 450*milli*darcy, 0.2);
rate   = .6*sum(poreVolume(G,rock))/T;
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
hT = computeTrans(G, rock);

%% Compute reference solution
xr = initState(G,W,100*barsa, [0 1]);
for n=1:128
    xr = incompTPFA(xr, G, hT, fluid, 'wells', W);
    xr = explicitTransport(xr, G, T/128, rock, fluid, 'wells', W);
end
pv = poreVolume(G,rock);
l1s = norm(xr.s(:,1).*pv,1);
l2p = norm(xr.pressure.^2.*pv, 1);

%% Contrast multiphase displacement with single-phase time lines
% For Corey exponent 2 and equal viscosities, the displacement front will
% propagate at a speed of a=1/(2(sqrt(2)-1)) relative to the total
% velocity. Compute the countour showing the position of the displacement
% front at time 0.6 PVI.
a  = 1/(2*(sqrt(2)-1));
x  = initState(G,W,100*barsa, [0 1]);
x  = incompTPFA(x, G, hT, fluid, 'wells', W);

tau = computeTimeOfFlight(x, G, rock,  'wells', W);

figure('Position',[300 550 1100 300]);
C = contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a, '-k','LineWidth',1); clf


%% Run simulation with different time steps
% Finally, we compute the multiphase flow solution at time 0.6 PVI using a
% sequence of increasing splitting time steps. As the number of splitting
% steps increases, the approximate solution gradually approaches the
% correct solution computable on this grid
nval = 10;
cval = linspace(0,1,nval+1); cval=.5*cval(1:end-1)+.5*cval(2:end);
colormap(flipud(.5*jet(nval)+.5*ones(nval,3)));
err = nan(7,2);
for n=0:6
    nstep = 2^n;
    x  = initState(G,W,100*barsa, [0 1]);
    for i=1:nstep
        x  = incompTPFA(x, G, hT, fluid, 'wells', W);
        x  = explicitTransport(x, G, T/nstep, rock, fluid, 'wells', W);
    end
    
    err(n+1,1) = norm((xr.s(:,1) - x.s(:,1)).*pv, 1)/l1s;
    err(n+1,2) = norm((xr.pressure - x.pressure).^2.*pv, 1)/l2p;

    if rem(n,2), continue, end
    
    subplot(1,4,n/2+1);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    hold on, plot(C(1,2:end),C(2,2:end),'k-','LineWidth',1); hold off
    axis equal; axis([0 domain(1) 0 domain(2)]); caxis([0 1]);
    title([num2str(nstep) ' steps'])
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
end

%% Show convergence table
convrate = log(bsxfun(@rdivide,err(1:end-1,:),err(2:end,:)))./log(2);
convrate = [nan nan; convrate];
clc
fprintf('\n\t  L1(S)     rate\t  L2(p)     rate\n');
fprintf('\t-----------------------------------------\n');
fprintf('\t%.3e   %.2f \t%.3e   %.2f\n', ...
    [err(:,1),convrate(:,1),err(:,2),convrate(:,2)]');
fprintf('\t-----------------------------------------\n');
