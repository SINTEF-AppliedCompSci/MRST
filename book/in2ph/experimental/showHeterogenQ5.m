%% Heterogeneous quarter five-spot
% In this example, we will study a heterogeneous quarter five-spot with
% petrophysical data generated from a random Gaussian field.
mrstModule add incomp

%% Set up model and parameters
% Use a fixed seed for the random generator so that results are
% reproducible from one run to the next.
load rndseed.mat; rng(S);
T        = 1.5;
dims     = [60 60];
domain   = [600 600]*meter;
G        = computeGeometry(cartGrid(dims,domain));
b = log(milli*darcy);
a = (log(darcy)-b)/(.4 - .05);
p = gaussianField(G.cartDims, [0.05 0.4], [3 11], 4.5);
K = exp(a*(p-.05)+b);
rock.poro = p(:);
rock.perm = K(:);
gravity reset off

rate   = sum(poreVolume(G,rock));
W = addWell([],G, rock, 1, 'Type', 'rate', ...
    'Val', rate, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
W = addWell(W,G, rock, G.cells.num, 'Type', 'rate', ...
    'Val', -rate, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);
hT = computeTrans(G, rock);
x0 = initState(G,W,100*barsa, [0 1]);

cval = linspace(0,1,11); cval=.5*cval(1:end-1)+.5*cval(2:end);

%% Compare fingering effects
% We start by comparing fingering effects in a heterogeneous quarter
% five-spot problem as function of viscosity ratio. We show that an
% unfavorable mobility ratio (viscosity of injected water is less than the
% viscosity of the resident oil) leads to pronounced fingers whereas a
% piston-like displacement with minimal fingering effects is observed for a
% favorable mobility ratio.
%
% To this end, we compute solutions for three different viscosity ratios,
% mu_w : mu_o = N, for N=10, 1, and 1/10. To compare the resulting
% displacement fronts, we use an estimate of 1D wave speeds to estimate how
% fast each displacement front propagates relative to a passive particle,
% and then use this to scale the final time so that the fronts have swept a
% comparable spatial region of the reservoir.
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
mb      = zeros(M,3);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    [x,ip] = deal(x0,1);
    for i=1:M
        x  = incompTPFA(x, G, hT, fluid, 'wells', W);
        x  = implicitTransport(x, G, dT, rock, fluid, 'wells', W);
        mb(i,n) = sum(x.s(:,2).*poreVolume(G,rock));
        dt(i) = dT;
        wellSol{i,n} = getWellSol(W,x, fluid);
        % ws = getWellSol(W,x, fluid);
        % for m=1:numel(ws), wellSol{i,n,m}=ws(m); end
        if ~mod(i,M/10) && ip<=5
            subplot(3,5,is); ip=ip+1; is=is+1;
            plotData(x); axis equal; axis([0 domain(1) 0 domain(2)]);
            set(gca,'XTick',[],'YTick',[]);
            if n==1,
                title(sprintf('t=%.2f PVI', sum(dt)));
            else
                set(gca,'position',get(gca,'position')+(n-1)*[0 .05 0 0]);
            end
            drawnow;
        end
    end
end

%% Plot well responses
ws = cellfun(@(x) x, wellSol, 'UniformOutput', false);
t  = cumsum(dt);
figure;
plot(t,cellfun(@(x) x(2).Sw, wellSol));   set(gca,'XLim',[0 1.5]);
figure;
plot(t,cellfun(@(x) x(2).wcut, wellSol)); set(gca,'XLim',[0 1.5]);
figure;
spv = sum(x0.s(:,2).*rock.poro.*G.cells.volumes);
plot(t,cumsum(bsxfun(@times, abs(cellfun(@(x) x(2).qOs, wellSol)), dt)));
hold on;
plot(t,spv-mb,'-.','LineWidth',3);
plot([0 T],[spv spv],'--k');
hold off;
set(gca,'XLim',[0 1.5]);