%% Inverted Gravity Column
% We consider a setup with a heavier fluid placed on top of a lighter
% fluid. The lighter fluid will move upward and the heavier fluid will move
% downward and gradually we will approach a stable steady state.
mrstModule add incomp

%% Define the model
gravity reset on
G     = computeGeometry(cartGrid([1, 1, 40], [1, 1, 10]));
rock  = makeRock(G, 0.1*darcy, 1);
fluid = initCoreyFluid('mu' , [0.30860, 0.056641]*centi*poise     , ...
                       'rho', [ 975.86,  686.54]*kilogram/meter^3, ...
                        'n' , [      2,       2],...
                        'sr', [     .1,      .2],...
                        'kwm',[  .2142,    .85]);

%% Initialize problem
% Set up the fluid distribution, compute pressure
hT = computeTrans(G, rock);
xr = initResSol(G, 100.0*barsa, 1.0); xr.s(end/2+1:end) = 0.0;
xr = incompTPFA(xr, G, hT, fluid); x0 = xr;

% Plot the 3D solution
figure(1);
subplot(1,3,1:2);
plotCellData(G,xr.s(:,1),'EdgeColor','none');
set(gca,'XTick',[],'YTick',[],'Ztick',[]); box on; view(3); caxis([0 1]);
title(sprintf('Time: %.2f years', 0));
colormap(flipud(parula));

% Plot 2D profile of the solution
subplot(1,3,3); 
plot(xr.s(:,1),G.cells.centroids(:,3),'-o',...
    'MarkerSize',8,'MarkerFaceColor',[.5,.5,.5]);
set(gca,'XTick',[],'YTick',[],'YDir','reverse');

%n=0; print('-depsc2',['inv-column-' num2str(n) '.eps']);

figure(2);
col = jet(16); k=1;
plot(G.cells.centroids(:,3), xr.pressure,'-o','MarkerSize',2,'color',col(k,:),'LineWidth',1.5);

%% Solve the problem
% We use a number of time steps to march the solution towards steady state.
% Since we do not have any movement in the lateral directions, we can use
% the explicit solver as a pure gravity segregation solver.
dt = 5*day; t=0; 
S = zeros(G.cells.num,151); S(:,1) = xr.s(:,1);
for i=1:150
    
    % Saturation step
    xr = explicitTransport(xr, G, dt, rock, fluid, 'onlygrav', true);
    t = t+dt;  
    
    % Plot 3D solution
    figure(1);
    subplot(1,3,1:2); cla
    plotCellData(G,xr.s(:,1),'EdgeColor','none'), box on; view(3);
    set(gca,'XTick',[],'YTick',[],'Ztick',[]); caxis([0 1]);
    title(sprintf('Time: %.2f years', t/year));
    
    % Plot 2D profile of the saturation
    subplot(1,3,3);
    plot(xr.s(:,1),G.cells.centroids(:,3),'-o',...
        'MarkerSize',8,'MarkerFaceColor',[.5,.5,.5]);
    set(gca,'XTick',[],'YTick',[],'YDir','reverse');
    drawnow; pause(0.1);
    
    %if any(i==[25 50 75 100 150])
    %  n=n+1; print('-depsc2',['inv-column-' num2str(n) '.eps']);  
    %end
    
    if ~mod(i,10)
        k = k+1;
        figure(2); hold on;
        plot(G.cells.centroids(:,3), xr.pressure,'-','Color',col(k,:),'LineWidth',1.5);
        hold off
    end
    
    % Compute pressure for the next step
    xr = incompTPFA(xr, G, hT, fluid);
    S(:,i+1)=xr.s(:,1);
end
figure(2);
hold on
plot(G.cells.centroids(:,3), xr.pressure,'+','Color',col(k,:),'MarkerSize',6);
hold off

%%
figure;
pcolor(S'); shading interp;
colormap(flipud(parula));
set(gca,'YDir','reverse'); 
box on; %set(gca,'XTick',[],'YTick',[]);