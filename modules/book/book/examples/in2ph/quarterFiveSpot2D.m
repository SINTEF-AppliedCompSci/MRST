%% Homogeneous quarter five-spot
% In this example, we compare and contrast incompressible single-phase and
% two-phase flow for a homogeneous quarter five-spot. Through the example,
% you will also be introduced to the way well solutions are represented in
% more advanced multiphase simulators based on the AD-OO framework.

mrstModule add incomp diagnostics

%% Set up model
% Square domain with homogeneous rock properties, no flow across the
% boundaries, no gravity, injection in the southeast and production in the
% northwest corners, both operating at fixed bottom-hole pressure
gravity reset off
cartDim = [128 128 1];
fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000,  850] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
domain = [250 250 20];
G      = computeGeometry(cartGrid(cartDim,domain));
rock   = makeRock(G, 100*milli*darcy, 0.2);
pv     = poreVolume(G, rock);

W = addWell([],G, rock, 1, 'Type', 'bhp', ...
    'Val', 100*barsa, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'bhp', ...
    'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
hT = computeTrans(G, rock);

% Figure, colormap and contour values
figure('Position',[300 550 1100 650]);
nval = 20;
cval = linspace(0,1,nval+1); cval=.5*cval(1:end-1)+.5*cval(2:end);
colormap(flipud(.5*jet(nval)+.5*ones(nval,3)));

%% Compute solution and plot saturation evolution
% To study the solution, we will plot saturation profiles at four instances
% in time up to slightly beyond water breakthrough; this corresponds
% dimensionless times 0.2 to 0.8 PVI. To get the solution provies as
% accurate as possible, we use the explicit transport solver and M substeps
% to advance the solution 0.2 PVI forward in time. We continue computing
% the solution up to time 1.2 PVI, but do not show snapshots of the
% saturation field for the two last time intervals.

% Compute an initial single-phase pressure solution, from which we estimate
% the final time that corresponds to 1.2 PVI if this flow field remains
% unchanged. With quadratic relperm curves and equal viscosities, the
% multiphase displacement front will propagate at a speed of
% a=1/(2(sqrt(2)-1)) relative to the total velocity.
x  = initState(G,W,100*barsa, [0 1]);
x  = incompTPFA(x, G, hT, fluid, 'wells', W);
T  = 1.2*sum(pv)/x.wellSol(1).flux;
a  = 1/(2*(sqrt(2)-1));

% Compute time-of-flight for the single-phase flow field and record the
% corresponding breakthrough time in the producer.
tau = computeTimeOfFlight(x, G, rock,  'wells', W);
tbf = tau(W(2).cells,1);

% Initialize number of time intervals, cell array to hold well solutions,
% and array to hold the oil in place
[N,M]    = deal(6,10);
wellSols = cell(N*M+1,1);  wellSols{1} = getWellSol(W, x, fluid);
oip      = zeros(N*M+1,1); oip(1) = sum(x.s(:,2).*pv);
for n=1:N
    fprintf(1,'Main step %d: ',n);
    for m=1:M
       x  = incompTPFA(x, G, hT, fluid, 'wells', W);
       x  = explicitTransport(x, G, T/(N*M), rock, fluid, 'wells', W);
       
       wellSols{(n-1)*M+m+1} = getWellSol(W, x, fluid);
       oip((n-1)*M+m+1) = sum(x.s(:,2).*pv);
       if x.s(W(2).cells,1)<eps
           W(2).bt = (n-1)*M+m+1;
       end
       fprintf(1,'%d, ',m);
    end
    fprintf(1,'\n');
    if n>4, continue, end
    
    % Plot multiphase solution
    subplot(2,4,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [0 cval 1], 'EdgeColor','none');
    hold on;
    
    % Plot corresponding time lines from single-phase solution
    contour(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(tau/T,G.cartDims), a*n/N, '-k','LineWidth',1);
    caxis([0 1]);
    axis equal; axis([0 domain(1) 0 domain(2)]);
    title(sprintf('t=%.2f PVI',n*.2));
    set(gca,'XTick',[],'YTick',[]);
    
    % Plot multiphase solution as function of single-phase time-of-flight
    subplot(2,4,4+n)
    set(gca,'position',get(gca,'position')+[0 .12 0 0]);
    plot(tau(:,1)/tbf,x.s(:,1),'.k','MarkerSize',4);
    set(gca,'XLim',[0 2]); drawnow;
end

%% Plot production curves
% First we plot the saturation in the perforated cell and the corresponding
% water cut, i.e., the fractional flow evaluated in the completion
dt = [0; ones(N*M,1)*T/(N*M)]; t = 1.2*cumsum(dt)/T;
figure;
plot(t,cellfun(@(x) x(2).Sw, wellSols),'--', ...
    t,cellfun(@(x) x(2).wcut, wellSols),'-','LineWidth',1);
legend('Sw in completion','Water cut','Location','NorthWest');
axis([0 max(t) -.05 1.0]);

%%
% Second, we plot the oil rate used in our simulation in units m^3/day
figure;
qOs = cellfun(@(x) abs(x(2).qOs), wellSols);
stairs(t, qOs([2:end end])*day,'LineWidth',1); axis([0 1.2 30 260]);

%%
% Last, we plot the cumulative oil production computed from the well
% solution and compare this with the amount of extracted oil derived from a
% mass-balance computation (initial oil in place minus current oil in
% place). We also include a horizontal line indicating the initial oil in
% place, and a straight line showing the amount of oil we would have
% extracted if oil was produced at the constant initial rate
figure
plot(t,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSols)), dt)));
hold on;
plot(t,oip(1)-oip,'-.','LineWidth',3);
plot([0 1.2],oip([1 1]),'-k',t,min(t*oip(1),oip(1)),'-k', ...
    t(W(2).bt+[0 0]),[0 oip(1)],'--k');
hold off; axis tight; axis([0 max(t) 0 1.05*oip(1)]);


%% Plotting with GUI from ad-core
mrstModule add ad-core
plotWellSols(wellSols,cumsum(dt))
