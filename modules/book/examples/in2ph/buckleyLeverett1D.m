%% Classic test case: 1D Buckley-Leverett
% We investigate the effect of time step on the accuracy and convergence of
% the implicit transport solver
mrstModule add incomp

G = computeGeometry(cartGrid([100,1]));
rock = makeRock(G, 100*milli*darcy, 0.2);

fluid = initSimpleFluid('mu' , [   1,    1] .* centi*poise     , ...
                        'rho', [1000, 1000] .* kilogram/meter^3, ...
                        'n'  , [   2,    2]);
bc = fluxside([], G, 'Left',   1, 'sat', [1 0]);
bc = fluxside(bc, G, 'Right', -1, 'sat', [0 1]);

hT = computeTrans(G, rock);
rSol = initState(G, [], 0, [0 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'bc', bc);

% Explicit tranport solver
rSole = explicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'verbose', true);

% Implicit transport solver: try with one time step
[rSoli, report] = ...
    implicitTransport(rSol, G, 10, rock, fluid, 'bc', bc, 'Verbose', true);

%% Load and display convergence history
% Here we have cheated a little: That is, we have copied the screendump to
% a file and used a text editor to manipulate it so that it can easily be
% reloaded
figure;
d = load('screendump.dat');
i = isnan(d(:,1)); d(:,1) = 1; d(i,1)=0; d(:,1)=cumsum(d(:,1));
plot(d(:,3),d(:,1),'o-','MarkerFaceColor',[0.5,0.5,0.5]);
set(gca,'XScale','log','YDir','reverse','XDir','normal'); axis tight
set(gca,'XTick',[1e-10 1e-5 1],'FontSize',12);
hold on;
j = find(i);
hold on
plot(repmat([8e-11 10],numel(j),1)',[d(j,1)+.5 d(j,1)+.5]','--r');
hold off
set(gcf,'Position',[1120 55 230 760],'PaperPositionMode','auto');

%% Plot results with various number of time steps
figure
plot(G.cells.centroids(:,1), rSole.s(:,1),'k--','LineWidth',1.5);
leg = cell(7,1); leg{1} = 'Expl: 199 steps';

n   = [4 10 20 40 100 200];
its = [0 0 0 0 0 0 0];
col = 'rgbcmk';
hold on
for k=1:numel(n)
    rSolt = rSol;
    for i=1:n(k)
        [rSolt, report] = ...
            implicitTransport(rSolt, G, 10/n(k), rock, fluid, 'bc', bc);
        its(k) = its(k) + report.iterations + report.vasted_iterations;
    end
    plot(G.cells.centroids(:,1),rSolt.s(:,1), [col(k) '-'],'LineWidth',1.5);
    leg{k+1} = sprintf('n=%3d: %3d its',n(k),its(k));
end
hold off
legend(leg{:});

