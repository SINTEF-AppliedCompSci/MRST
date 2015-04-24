%% Gradient-based search for optimium of a quadratic function on unit square

% Define test-function parameters:
opts = {'th', pi/10, ...
        'rx', 4 ,    ...
        'ry', 1 ,    ...
        'x0', .8,    ...
        'y0', .8};
% Compute function value at mesh-points for subsequent contour-plots
[X,Y] = meshgrid(0:.01:1, 0:.01:1);
Z = testFunction(X,Y, opts{:});
% Function handle:
f = @(u)testFunction(u(1),u(2), opts{:});
% Initial guess:
u0 = [.1, .1]';

% Linear inequality constaint in addition to box: y< -x + 1.2
le = struct('A', [1 1], 'b', 1.2);
%% Optimize using BFGS:
[v1,u1,hst1] = unitBoxBFGS(u0, f, 'useBFGS', true, ...
                                   'wolfe1', 1e-3, ...
                                   'wolfe2',  0.5,  ...
                                   'linIneq',  le);


% Gather evolution of control-vector:
p1 = horzcat(hst1.u{:})';
% Plot contour of objective and evolution of controls
figure(11),hold on
contour(X,Y,Z,-.5:.011:0);
plot(p1(:,1),p1(:,2),'.-g','LineWidth', 2, 'MarkerSize', 20)
plot(p1(end,1),p1(end,2),'or','LineWidth', 2, 'MarkerSize', 15)
% Plot linear inequality constraint
x = 0:1; y = -(le.A(1)*x - le.b)/le.A(2);
plot(x, y, '--k','LineWidth', 2)
axis equal, axis tight, box on
axis equal, axis tight, box on, axis([0 1 0 1])
