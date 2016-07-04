%% Quarter five-spot: Illustration of temporal splitting errors
% We use the standard quarter five-spot test case to illustrate splitting
% errors that may arise when using a sequential solution procedure to
% simulate multiphase flow 
mrstModule add incomp

%% Set up model
T = 5*year;
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

cval = linspace(0,1,11); cval=.5*cval(1:end-1)+.5*cval(2:end);

figure('Position',[300 550 1100 300]);
colormap(flipud(.5*jet(10)+.5*ones(10,3)));

%% Plot of evolving displacement front
% First, we use contour lines to show four snapshots of the evolving
% displacement front. To get the front as accurately as possible, we use an
% explicit scheme and update the pressure for (almost) every step of the
% transport solver.

% NB: REMEMBER TO REPORT THE NUMBER OF TIME STEPS!

x  = initState(G,W,100*barsa, [0 1]);
for n=1:4
    for i=1:8
       x  = incompTPFA(x, G, hT, fluid, 'wells', W);
       x  = explicitTransport(x, G, T/24, rock, fluid, 'wells', W);
    end
    subplot(1,4,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    axis equal; axis([0 domain(1) 0 domain(2)]);
    title(sprintf('t=%.2f PVI',n*.2));
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
end

%% Contrast multiphase displacement with single-phase time lines
% For Corey exponent 2 and equal viscosities, the displacement front will
% propagate at a speed of a=1/(2(sqrt(2)-1)) relative to the total
% velocity. To demonstrate the difference between single-phase and
% multiphase, we therefore add time lines at time-of-flight equal
% a*t. 
a  = 1/(2*(sqrt(2)-1));
x  = initState(G,W,100*barsa, [0 1]);
x  = incompTPFA(x, G, hT, fluid, 'wells', W);
mrstModule add diagnostics
tau = computeTimeOfFlight(x, G, rock,  'wells', W);
for i=1:4
    subplot(1,4,i);
    hold on;
    contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a*i/3, '-k','LineWidth',1);
    caxis([0 1]);
    hold off
end

%% Run simulation with different time steps
% Finally, we consider the solution at time 0.6 PVI and compute splitting
% solutions with a sequence of increasing splitting time steps. As the
% number of splitting steps increases, the approximate solution gradually
% approaches the correct solution.
nstep = [1 4 16 64];
figure('Position',[300 550 1100 300]);
C= contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a, '-k','LineWidth',1);
colormap(flipud(.5*jet(10)+.5*ones(10,3)));
for n=1:numel(nstep)
    x  = initState(G,W,100*barsa, [0 1]);
    for i=1:nstep(n)
        x  = incompTPFA(x, G, hT, fluid, 'wells', W);
        x  = explicitTransport(x, G, T/nstep(n), rock, fluid, 'wells', W);
    end
    subplot(1,4,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    hold on, plot(C(1,2:end),C(2,2:end),'k-','LineWidth',1); hold off
    axis equal; axis([0 domain(1) 0 domain(2)]); caxis([0 1]);
    title([num2str(nstep(n)) ' steps'])
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
end