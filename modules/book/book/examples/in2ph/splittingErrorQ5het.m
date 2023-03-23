%% Heterogeneous quarter five-spot: temporal splitting errors
% In this example, we use a heterogeneous quarter five-spot with
% petrophysical data sampled from the topmost layer of the SPE 10 data set
% as an example to illustrate splitting errors that may arise when using a
% sequential solution procedure. We consider three different mobility
% ratios: an unfavorable, equal mobilities, and favorable.

mrstModule add incomp spe10 ad-core

%% Set up model and parameters
dims      = [60 120 1];
domain    = dims.*[20 10 2]*ft;
G         = computeGeometry(cartGrid(dims,domain));
rock      = getSPE10rock((1:dims(1)),(1:dims(2))+60,1);
rock.poro = max(rock.poro,.0005);
pv        = poreVolume(G,rock);
gravity reset off

rate   = sum(poreVolume(G,rock));
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
hT = computeTrans(G, rock);
x0 = initState(G,W,100*barsa, [0 1]);


%% Self convergence with respect to time step
% For the first test, we illustrate difference in self convergence as a
% function of mobility ratio M = muw/muo for M=0.1, 1, and 10. These cases
% have a leading displacement profile that moves faster than the Darcy
% velocity by a factor rvel=1/(2M(sqrt((M+1)/M)-1)). For each mobility
% ratio, we set up a simulation that would take 2^m time splitting steps to
% reach 0.5 PVI for m=1:6. However, to get saturation profiles that cover a
% more similar portion of the reservoir, we scale the end time for M=0.1
% and M=10 by rvel(1)/rvel(M). The self convergence is measured relative to
% a solution that would have used 256 steps to get to 0.5 PVI.
mu     = [1 10; 1 1 ; 10 1];
rvel   = @(M) 1./(2.*M.*(sqrt((M+1)./M)-1));
tscale = rvel(mu(:,1)./mu(:,2));
tscale = tscale(2)./tscale;
T      = 0.5;

% Prepare plotting
set(gcf,'Position',[200 550 1240 480]);
cval = linspace(0,1,11); 
cval = .5*cval(1:end-1)+.5*cval(2:end);
plotData = @(x) ... 
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
    reshape(G.cells.centroids(:,2), G.cartDims), ...
    reshape(x.s(:,1),G.cartDims), [cval 1]);
colormap(flipud(.5*jet(10)+.5*ones(10,3)));

% Loop over all three viscosity ratios
err = zeros(3,6);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    for m=[8 1:6]
        
        % Time loop
        [x,dt,t,tend]  = deal(x0, T*2^(-m), 0, tscale(n)*T);
        while t<tend
            x = incompTPFA(x, G, hT, fluid, 'wells', W);
            [tl,dtl,tle] = deal(0,T/2^8, min(tend-t,dt));
            while tl < tle
                x = implicitTransport(x, G, dtl, rock, fluid, 'wells', W);
                tl = tl+dtl;
            end
            t = t+dt;
        end
        
        % Compute discrepancy from reference solution
        if m==8
            sref = x.s(:,1); snorm = sum(abs(sref).*pv);
        else
            err(n,m) = sum(abs(x.s(:,1)-sref).*pv)/snorm;
        end
        
        % Plot contour plot of solution
        subplot(3,7,sub2ind([7,3],min(m,7),n))
        plotData(x); axis equal; axis([0 domain(1) 0 domain(2)]);
        if n==1
            title(sprintf('%d steps',2^m));
        else
            pos=get(gca,'position');
            set(gca,'position',pos+(n-1).*[0 .05 0 0]);
        end
        set(gca,'XTick',[],'YTick',[]);
        drawnow;
    end
end

% Plot convergence
figure
semilogy(err','o-','MarkerSize',7,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
set(gca,'XTick',1:6,'XTickLabel',{2.^(1:6)}); legend('1:10','1:1','10:1');

%% Well curves convergence
% We repeate the same simulation as above, except that we continue until
% 1.5 PVI for all three mobility ratios. We run simulations with constant
% time-step length (3*4^[0:4] steps) and with a rampup.
Tn    = 1.5;

% Constant time step
nstep = [1 4 16 64 256];
dt    = T/nstep(end);
wellSols = cell(nstep(end),numel(nstep),3);

% Rampup 
dT = repmat(T/32, Tn/T*32,1);
dT = [dT(1)*sort(repmat((2.^-[1:4 4])',2,1)); dT(3:end)];
wSol     = cell(numel(dT),3);
for n=1:3
    fprintf('Fluid %d: ', n);
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    for k=1:numel(nstep)
        x = x0;
        ws = 1;
        fprintf('%d,', nstep(k));
        for i=1:nstep(k)*Tn/T
            x = incompTPFA(x, G, hT, fluid, 'wells', W);
            for m=1:nstep(end)/nstep(k)
                x = implicitTransport(x, G, dt, rock, fluid, 'wells', W);
                wellSols{ws,k,n} = getWellSol(W, x, fluid); ws = ws+1;
            end
        end
    end
    fprintf('rampup\n');
    x = x0;
    for i=1:numel(dT)
        x = incompTPFA(x, G, hT, fluid, 'wells', W);
        x = implicitTransport(x, G, dT(i), rock, fluid, 'wells', W);
        wSol{i,n} = getWellSol(W, x, fluid);
    end
end


%%
n = 2;
[tws,dsn,t] = deal(cell(1,numel(nstep)+1));
for i=1:numel(nstep)
    tws(i) = {vertcat(wellSols(:,i,n))};
    t(i)   = {cumsum(dt*ones(Tn/T*nstep(end),1))};
    dsn(i) = {num2str(nstep(i))};
end
tws(end) = {vertcat(wSol(:,n))};
t(end)   = {cumsum(dT)};
dsn(end) = {'rampup'};

plotWellSols(tws,t, 'datasetnames', dsn);