%% Section 5: Adapting cell centroids
% This section discusses how to adapt grids to cell-centroid constraints.
colors = get(gca,'colororder');

%% Adaption with/without protection layer
% The first example compares adaption to a piecewise linear path with and
% without the use of protection sites, which are introduced so that the
% adapted cells become more rectangular. (See Figure 10)
cellConstraints = {[0, 0.4; 0.2, 0.5; 0.6, 0.5; 0.8, 0.6]};
[X, Y]   = meshgrid(linspace(0,1,10));
distance = @(x) 0.12 * ones(size(x, 1), 1);

for i=1:2
    % Add background sites to the domain
    [CCsites, cGs, protSites, pGs] = lineSites2D(cellConstraints, 0.12, ...
        'protLayer',logical(1-i),'protD', {distance});
    bgSites = removeConflictPoints([X(:), Y(:)], CCsites, cGs);
    bgSites = removeConflictPoints(bgSites, protSites, pGs);
    sites   = [CCsites; protSites; bgSites];
    bnd     = [0, 0; 1, 0; 1, 1; 0, 1];
    G       = clippedPebi2D(sites, bnd);

    % Plot result
    figure(i),clf
    G = computeGeometry(G);
    plotGrid(G, 'facecolor', 'none')
    hold on
    plotLinePath(cellConstraints,'-','color',colors(4, :), 'linewidth', 3);
    plot(CCsites(:, 1), CCsites(:,2), '.','markersize',35,'Color',[.5 .5 .7])
    plot(CCsites(:, 1), CCsites(:,2), '.k','markersize',10)
    plot(G.cells.centroids(1:size(CCsites, 1), 1), ...
        G.cells.centroids(1:size(CCsites, 1), 2),...
        'xr', 'markersize',10,'linewidth', 2)
    plot(protSites(:, 1), protSites(:, 2),'.','color',colors(2,:),'markersize',15)
    axis equal tight off
    hold off
end

%% Adapt to cell and face centroids
% The next example adds wells to the example from Section 4.1 so that the
% grid adapts to constrains on both faces and cell centroids (as shown in
% Figure 12 in the chapter)
lines = {[0.2, 0.2; 0.7, 0.05], ...
         [0.2, 0.05; 0.7, 0.2], ...
         [0.1, 0.4; 0.6, 0.6], ...
         [0.1, 0.7; 0.45, 0.7; 0.55, 0.3]};
wellLines = {[0.1, 0.6; 0.2, 0.6; 0.3, 0.5; 0.4, 0.3], [0.8, 0.8]};

G = pebiGrid2D(0.06, [1, 1], 'faceConstraints', lines, ...
    'cellConstraints', wellLines, ...
    'CCRefinement',true, 'CCFactor', 0.25);

% Plot the resulting grid
figure
subplot(1,3,2:3)
plotGrid(G,'facecolor','none');
color = get(gca, 'colororder');
plotGrid(G, G.cells.tag, 'facecolor', color(3,:))
plotFaces(G, find(G.faces.tag),'edgecolor', color(2, :),'LineWidth',2)
plotLinePath(lines, ':o', 'color', color(1, :), 'LineWidth',4,'MarkerSize',3);
plotLinePath(wellLines(1), ':o', 'color', color(4, :), 'LineWidth', 4,'MarkerSize',3)
axis equal tight off
ax1 = gca;

subplot(2,3,4)
copyobj(get(ax1,'Children'),gca);
axis([.25 .35 .43 .48]);
ax=gca;
ax.Position(2)=0.177;
ax.Position(1) = 0.17;
axis equal, box on, set(gca,'XTick',[],'YTick',[]);

%% Controlling the adaption
% The last example shows a few parameters you can use to control how the
% grid adapts to cell constraints
w     = {[0.2,0.3;0.5,0.5;0.8,0.5]};
gS    = [0.1,0.1];
pdims =[1,1];

figure
col = get(gca,'colororder');
for i=1:4
    switch i
        case 1
            pargs = {'CCFactor',1/2,'mlqtMaxLevel',1};
            name  = 'mlqtMaxLevel=1';
        case 2
            pargs = {'CCFactor',1/16,'mlqtMaxLevel',3};
            name  = 'mlqtMaxLevel=3';
        case 3
            pargs = {'CCFactor',1/4, 'mlqtMaxLevel',2,'mlqtLevelSteps',[.05,.1]};
            name  = 'mlqtLevelSteps=[.05,.1]';
        case 4
            pargs = {'CCFactor',1/4, 'mlqtMaxLevel',2,'mlqtLevelSteps',[.1,.2]};
            name  = 'mlqtLevelSteps=[.1,.2]';
    end

    subplot(2,2,i),cla
    G = compositePebiGrid2D(gS, pdims, 'cellConstraints', w, pargs{:});
    plotGrid(G,'FaceColor','none'); axis equal tight off
    plotGrid(G,G.cells.tag,'facecolor',col(3,:))
    plotLinePath(w,'--o','linewidth',2,'color',col(4,:),'MarkerFaceColor','w','MarkerSize',4);
    title(name,'FontWeight','normal');
end

