%% Heterogeneous quarter five-spot
% In this example, we will study a heterogeneous quarter five-spot with
% petrophysical data sampled from the topmost layer of the SPE 10 data set.

mrstModule add incomp spe10

%% Set up model and parameters
T        = 0.5;
dims     = [60 120 1];
domain   = dims.*[20 10 2]*ft;
G        = computeGeometry(cartGrid(dims,domain));
rock     = SPE10_rock((1:dims(1)),(1:dims(2))+60,1);
rock.perm = rock.perm*milli*darcy;
rock.poro = max(rock.poro,.0005);
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
N      = 10;
mu     = [1 N; 1 1 ; N 1];
g      = @(N) N./(2*(sqrt(N+1)-1));
tscale = g([N 1 1/N]);
tscale = tscale(2)./tscale;
M      = 40;
set(gcf,'Position',[300 550 1100 300]);
colormap(flipud(.5*jet(10)+.5*ones(10,3)));
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    x  = x0;
    for i=1:round(tscale(n)*M)
       x  = incompTPFA(x, G, hT, fluid, 'wells', W);
       x  = implicitTransport(x, G, T/M, rock, fluid, 'wells', W);
    end
    subplot(1,3,n);
    contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
        reshape(G.cells.centroids(:,2), G.cartDims), ...
        reshape(x.s(:,1),G.cartDims), [cval 1]);
    axis equal; axis([0 domain(1) 0 domain(2)]);
    title(sprintf('Viscosity ratio %d:%d',mu(n,:)));
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
end

%% Show evolution of 

%% Self-convergence with respect to time step
% We repeat the experiment above and measure the self-convergence with
% respect to number of time steps for the three different viscosity ratios
pv = poreVolume(G,rock);
err = zeros(3,6);
set(gcf,'Position',[300 550 1100 480]);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    for m=[7 1:6]
        x  = x0;
        M = round(tscale(n)*2^m);
        for i=1:M;
            x = incompTPFA(x, G, hT, fluid, 'wells', W);
            x = implicitTransport(x, G, tscale(n)*T/M, rock, fluid, 'wells', W);
        end
        if m==7,
            sref = x.s(:,1);
            snorm = sum(abs(sref).*pv);
        else
            err(n,m) = sum(abs(x.s(:,1)-sref).*pv)/snorm;
            subplot(3,6,sub2ind([6,3],m,n))
            contourf(reshape(G.cells.centroids(:,1), G.cartDims),...
               reshape(G.cells.centroids(:,2), G.cartDims), ...
               reshape(x.s(:,1),G.cartDims), [cval 1]);
            axis equal; axis([0 domain(1) 0 domain(2)]);
            if n==1, 
                title(sprintf('%d steps',2^m));
            else
                pos=get(gca,'position');
                set(gca,'position',pos+(n-1).*[0 .025 0 0]);
            end
            set(gca,'XTick',[],'YTick',[]);
            drawnow;
        end
    end
end
figure
semilogy(err','o-','MarkerSize',7,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
set(gca,'XTick',1:6,'XTickLabel',{2.^[1:6]}); legend('1:10','1:1','10:1');

%% Well curves: error as function of time step
T = 1.5;
[Mx,My,Mz] = deal(2^3*T/0.5, 2^5*T/0.5, 2^7*T/0.5);
[tx,ty,tz] = deal(zeros(Mx,1),zeros(My,1),zeros(Mz,1));
[aerrx,rerrx,wx] = deal(zeros(Mx,3));
[aerry,rerry,wy] = deal(zeros(My,3));
wz = zeros(Mz,3);
spv  = sum(pv);
for n=1:3
    fluid = initSimpleFluid('mu', mu(n,:).* centi*poise, ...
        'rho', [1000,1000].*kilogram/meter^3, 'n',  [2 2]);
    [x,y,z] = deal(x0);
    [ij,ik] = deal(1);
    for i=1:Mx
        disp(i)
        x = incompTPFA(x, G, hT, fluid, 'wells', W);
        x = implicitTransport(x, G, T/Mx, rock, fluid, 'wells', W);
        wx(i,n) = x.s(W(2).cells,1);
        tx(i) = T/Mx;
        for j=1:My/Mx
            y = incompTPFA(y, G, hT, fluid, 'wells', W);
            y = implicitTransport(y, G, T/My, rock, fluid, 'wells', W);
            wy(ij,n) = y.s(W(2).cells,1);
            ty(ij) = T/My;
            for k=1:Mz/My
                z = incompTPFA(z, G, hT, fluid, 'wells', W);
                z = implicitTransport(z, G, T/Mz, rock, fluid, 'wells', W);
                wz(ik,n) = z.s(W(2).cells,1);
                tz(ik) = T/Mz;
                ik = ik+1;
            end
            aerry(ij,n) = sum(abs(z.s(:,1)-y.s(:,1)).*pv);
            rerry(ij,n) = aerry(ij,n)/sum(z.s(:,1).*pv);
            aerry(ij,n) = aerry(ij,n)/spv;
            ij = ij+1;
        end
        aerrx(i,n) = sum(abs(z.s(:,1)-x.s(:,1)).*pv);
        rerrx(i,n) = aerrx(i,n)/sum(z.s(:,1).*pv);
        aerrx(i,n) = aerrx(i,n)/spv;
    end
end
%%
figure
semilogy(cumsum(tx),aerrx,'o-','MarkerSize',7,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
hold on
semilogy(cumsum(ty),aerry,'-s','MarkerSize',4,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
hold off
legend('1:10, n=8','1:1','10:1','1:10, n=32','1:1','10:1'); axis tight
%
figure
semilogy(cumsum(tx),rerrx,'o-','MarkerSize',7,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
hold on
semilogy(cumsum(ty),rerry,'-s','MarkerSize',4,'MarkerFaceColor',[.8 .8 .8],'LineWidth',1);
hold off
legend('1:10, n=8','1:1','10:1','1:10, n=32','1:1','10:1'); axis tight
%
figure; 
plot(cumsum(tx),wx,'--o',cumsum(ty),wy,'-s',cumsum(tz),wz,'-','MarkerSize',4);