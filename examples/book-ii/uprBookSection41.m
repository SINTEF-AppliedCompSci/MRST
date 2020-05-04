%% Section 4.1: Simplex conformity to face constraints
% This script presents two examples that illustrate simplex conformity to
% face constraints. The first example imposes the constraints manually
% using low-level routines from the UPR module to demonstrate the basic
% principles. The second example shows how to do it automatically

%% Making a constrained grid manually
% Curve constraint with vertices
figure
colors = get(gca,'colororder');
W = {[0, 0.4; 0.2, 0.5],[0.2, 0.5; 0.4, 0.5],[0.4, 0.5; 0.6, 0.6]};  % Vertices
v = [0, 0.4; 0.2, 0.5; 0.4, 0.5; 0.6, 0.6];  % Vertices
plot(v(:,1),v(:,2),'-o','LineWidth',3,'Color',colors(1,:),'MarkerSize',5);
axis equal, axis([0 1 0 1])

%%
% Circle at each vertex
d = sqrt(sum(diff(v).^2,2));       % Distance between vertices
R = .6*min(d);                     % Circle radius
theta = linspace(0,2*pi)';
hold on
for j = 1:size(v,1)
    X = bsxfun(@plus, v(j,:), R*[cos(theta), sin(theta)]);
    plot(X(:,1), X(:,2),'k:', 'linewidth', 2)
end

%%
% Calculate circle intersections
dn = sqrt(R^2 - (d/2).^2);         % Normal offset
t = bsxfun(@rdivide, diff(v), d);  % Tangent vector
n = [-t(:, 2), t(:, 1)];           % Normal vector

% Add sites at both intersections
center      = v(1:end-1, :) + bsxfun(@times, d/2, t);
left_sites  = center + bsxfun(@times, dn, n);
right_sites = center - bsxfun(@times, dn, n);
sites = [left_sites; right_sites];
plot(sites(:, 1), sites(:, 2), '.', 'markersize', 25, 'color', colors(2,:))

% Add a site on the last circle
tip_sites = v(end, :) + R / sqrt(2);
plot(tip_sites(1), tip_sites(2), '.','markersize', 25, 'color', colors(5,:))
sites = [sites; tip_sites];

%%
% Add background sites to the domain
[X, Y]    = meshgrid(0:.2:1, 0:.2:1);
res_sites = removeConflictPoints([X(:), Y(:)], [v; tip_sites], R);
plot(res_sites(:, 1), res_sites(:, 2), '.', 'markersize', 25, 'color', [.6 .6 .6])
sites     = [sites; res_sites];

%%
% Create and plot the PEBI grid
bnd = [0, 0; 1, 0; 1, 1; 0, 1];
G = clippedPebi2D(sites, bnd);
plotGrid(G, 'facecolor', 'none')
hold off, axis off

%% Making a constrained grid using high-level library routines
% Define lower-dimensional lines and construct a composite grid with
% structured background reservoir sites
lines = {[0.2, 0.2; 0.7, 0.05], ...
         [0.2, 0.05; 0.7, 0.2], ...
         [0.1, 0.4; 0.6, 0.6], ...
         [0.1, 0.7; 0.45, 0.7; 0.55, 0.3]};
G = compositePebiGrid2D([0.05 .05], [1, 1], 'faceConstraints', lines);

% Uncomment the next line to get a grid with unstructured background sites
% instead:
% G = pebiGrid2D(.05,[1 1],'faceConstraints', lines);

%%
% Plot grid and show two zooms on specific details
figure('Position',[840 360 640 420])
subplot(1,3,2:3)
plotGrid(G,'facecolor','none');
color = get(gca, 'colororder');
plotFaces(G, find(G.faces.tag),'edgecolor', color(2, :),'LineWidth',2)
plotLinePath(lines,':o','color',color(1,:), 'LineWidth', 4, 'MarkerSize',3);
axis equal tight off
ax1 = gca;

subplot(2,3,4)
copyobj(get(ax1,'Children'),gca);
axis equal tight, axis([.35 .55 .03 .22]);
ax=gca;
ax.Position(2)=0.122;
ax.Position(1) = 0.17;
box on, set(gca,'XTick',[],'YTick',[]);

subplot(2,3,1), hold on
copyobj(get(ax1,'Children'),gca);
axis equal tight, axis([.1 .3 0. .16]);
ax=gca;
ax.Position(2)=0.595;
ax.Position(1) = 0.17;
plot(.1731,.0758,'sr','MarkerSize',24,'LineWidth',1)
box on, set(gca,'XTick',[],'YTick',[]);