%% Section 4.2: Configuring the simplex-conformity methods
% This section discusses how to configure the simple-conformity method
% introduced in uprBookSection41.m. The first two examples demonstrate how
% your can use interpolation of face constraints. The last example 

%% Interpolation versus exact representation of face constraints
% The first example has three constraints: a curved constraint represented
% by 100 points and two L-shaped constraints that each is represented by
% two line segments. The first two constraint are set to be interpolated,
% whereas the second L-shaped constraint is imposed exactly as given.
x     = linspace(0.2,0.8);
y     = 0.4+0.3*sin(pi*x);
lines = {[x(:), y(:)], ...
         [.0,.65; .0,.4; .25,.4],...
         [.75,.4; 1.0, .4; 1.0, .65]};
F = surfaceSites2D(lines, 0.1,'interpolateFC', [true, true, false]);

% Exercise: try to change the interpolation method for the curved
% constraint and see how this affects the grid.

%%
% Plot the tessellation of the constraints (in red color), the
% corresponding edge sites (orange) and the tip sites (purple).

subplot(2,1,1),cla, hold all
colors = get(gca,'colororder');
plotLinePath(lines,'--','LineWidth',2,'color',colors(1,:));
plot(F.c.CC(:,1), F.c.CC(:,2), '.','markersize',24,'color',colors(2,:))
plot(F.f.pts(:,1),F.f.pts(:,2),'.','markersize',24,'color',colors(3,:))
plot(F.t.pts(:,1),F.t.pts(:,2),'.','markersize',24,'color',colors(4,:))
axis off equal

% If we had not used interpolation, the tessellation of the curved
% constraints would have been as dense as the defining points. With
% interpolation, we only get eight points that are approximately
% equidistant along the curve. The drawback with interpolation is that the
% constraint will not be represented exactly. This is evident at the kink
% in the left L-shaped constraint, which is not reproduced by any of the
% faces in the resulting grid.

subplot(2,1,2),cla, hold all
[X, Y]    = meshgrid(-1:.1:2, -1:.1:2);
res_sites = removeConflictPoints([X(:), Y(:)], [F.f.pts;F.t.pts], ...
                                 [F.f.Gs; 0.04;0.04;0.04;0.04;0.04;0.04]);
sites     = [F.f.pts; F.t.pts; res_sites];
bnd = [-1, -1; 2, -1; 2, 2; -1, 2];
G = clippedPebi2D(sites, bnd);
plotLinePath(lines,':','color',colors(1,:), 'linewidth',3)
plotGrid(G, 'facecolor', 'none')
hold off, axis equal off
axis([-0.1, 1.1, 0.32, .78])

%% Adapting the density in the tessellation of the constraint
% The interpolation function also accepts a function handle that enables
% you tocontrol the discretization parameter in space. This is illustrated
% by the following example. To produce the plot in the upper subfigure, you
% must first enable plotting in the interLinePath function by commenting in
% again all lines marked by a comment '% makePlot'
figure
subplot(2,1,1)
F = surfaceSites2D(lines(1),0.1,'interpolateFC', true,...
    'distFun', @(x) .05+.125*x(:,1));
axis off, axis([.1 1 .45 .7])

subplot(2,1,2), hold all
colors = get(gca,'colororder');
plotLinePath(lines(1),'--','LineWidth',2,'color',colors(1,:));
plot(F.c.CC(:,1), F.c.CC(:,2), '.','markersize',24,'color',colors(2,:))
plot(F.f.pts(:,1),F.f.pts(:,2),'.','markersize',24,'color',colors(3,:))
plot(F.t.pts(:,1),F.t.pts(:,2),'.','markersize',24,'color',colors(4,:))
axis off equal

%% Controlling how the volumetric grid adapts to the fault constraint
% We compare six different setups to illustrate the different parameters
% you can use to control how the volumetric grid adapts to the face
% constraints. (This produces the plot in Figure 7 of the upr chapter). To
% avoid duplicating code, we use a local function for plotting, which will
% not work if you simply copy the code or execute it single line from
% outside of the script.
figure
f     = {[0.2,0.3;0.5,0.5;0.8,0.5]};

subplot(2,3,1),cla
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f));
title('Default setup','FontWeight','normal');

subplot(2,3,2), cla
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f, 'FCFactor',1/4));
title('FCFactor=1/4','FontWeight','normal');

subplot(2,3,3), cla
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f, 'circleFactor',0.95));
title('circleFactor=0.95','FontWeight','normal');

subplot(2,3,4)
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f, ...
                          'FCFactor',1/4,'FCRefinement',true));
title('FCRefinement','FontWeight','normal');

subplot(2,3,5)
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f, 'FCFactor',1/4, ...
    'FCRefinement',true, 'FCEps',.5));
title('FCEps=1/2','FontWeight','normal');

subplot(2,3,6)
myPlotGrid(f, pebiGrid2D(.1, [1 1],'faceConstraints',f, ...
    'FCRho', @(p) 1 - p(:,1)));
title('FCRho=@(p) 1 - p(:,1)','FontWeight','normal');

function myPlotGrid(f, G)
    plotGrid(G,'FaceColor','none'); axis tight off
    plotLinePath(f,'--o','color',[0 .447 .741], ...
        'linewidth',2,'MarkerFaceColor','w');
    axis equal tight
end