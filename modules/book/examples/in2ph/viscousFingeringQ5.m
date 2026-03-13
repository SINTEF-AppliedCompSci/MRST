%% Heterogeneous quarter five-spot
% In this example, we will study a heterogeneous quarter five-spot with
% petrophysical data sampled from the SPE 10 benchmark. We compare
% fingering effects in a quarter five-spot problem as function of viscosity
% ratio. We show that an unfavorable mobility ratio (viscosity of injected
% water is less than the viscosity of the resident oil) leads to pronounced
% fingers, whereas a piston-like displacement with minimal fingering effects
% is observed for a favorable mobility ratio.
mrstModule add incomp ad-core spe10

%% Set up model and parameters
T         = 20*year;
dims      = [60 120 1];
domain    = dims.*[20 10 2]*ft;
G         = computeGeometry(cartGrid(dims,domain));
rock      = getSPE10rock((1:dims(1)),(1:dims(2))+60,1);
rock.poro = max(rock.poro,.0005);
gravity reset off

rate   = sum(poreVolume(G,rock))/T;
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
hT = computeTrans(G, rock);
x0 = initState(G,W,100*barsa, [0 1]);

cval = linspace(0,1,11); cval=.5*cval(1:end-1)+.5*cval(2:end);

%% Compare fingering effects
% We compute solutions for three different viscosity ratios, mu_w : mu_o =
% N, for N=10, 1, and 1/10.
set(gcf,'Position',[250 150 1100 660]);
plotData = @(x) ... 
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
    reshape(G.cells.centroids(:,2), G.cartDims), ...
    reshape(x.s(:,1),G.cartDims), [cval 1]);
colormap(flipud(.5*jet(10)+.5*ones(10,3)));

[N,M]   = deal(10,200);
mu      = [1 N; 1 1 ; N 1];
[dt,dT] = deal(zeros(M,1), T/M);
is      = 1;
wellSol = cell(M,3);
oip     = zeros(M,3);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    [x,ip] = deal(x0,1);
    for i=1:M
        x  = incompTPFA(x, G, hT, fluid, 'wells', W);
        x  = implicitTransport(x, G, dT, rock, fluid, 'wells', W);
        
        dt(i) = dT;
        oip(i,n) = sum(x.s(:,2).*poreVolume(G,rock));
        wellSol{i,n} = getWellSol(W,x, fluid);

        if ~mod(i,40) && ip<=5
            subplot(3,5,is); ip=ip+1; is=is+1;
            plotData(x); axis equal; axis([0 domain(1) 0 domain(2)]);
            set(gca,'XTick',[],'YTick',[]);
            if n==1
                title(sprintf('%.0f years', sum(dt)/year));
            else
                set(gca,'position',get(gca,'position')+(n-1)*[0 .05 0 0]);
            end
            drawnow;
        end
    end
end

%% Plot well responses
% This plotting only works for predefined quantities like surface rates,
% bottom-hole pressure, etc. If the well structure contains other
% quantities, these will appear in the list of available fields but cannot
% be plotted.
plotWellSols({wellSol(:,1),wellSol(:,2),wellSol(:,3)},...
    cumsum(dt),'datasetnames',{'Ratio 1:10','Ratio 1:1','Ratio 10:1'});