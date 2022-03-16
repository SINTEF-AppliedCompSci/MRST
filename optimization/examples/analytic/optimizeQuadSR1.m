%% Unconstrained, gradient-based search for optimium of a quadratic function

mrstModule add ad-core optimization
% Define test-function parameters:
opts = {'th', pi/10, ...
        'rx', 4 ,    ...
        'ry', 1 ,    ...
        'x0', 2.4,    ...
        'y0', 0.3, ...
        'invertSign', true};

% Compute function value at mesh-points for subsequent contour-plots
x0 = opts{8}; y0 = opts{10};
[X,Y] = meshgrid(-1:x0/100:x0, -1:y0/100:y0);
Z = testFunction(X,Y, opts{:});

% Function handle:
f = @(u)testFunction(u(1),u(2), opts{:});

% Initial guess:
u0 = [.1, .1]';

%% Minimize using LR1:
[v, u, hst1] = optimizeSR1(u0, f, ...
                           'delta', 2, ...
                           'plotEvolution', true);


% Gather evolution of control-vector:
p1 = horzcat(hst1.u{:})';
% Plot contour of objective and evolution of controls
figure(11),hold on
contour(X,Y,Z,0:f(u0)/100:f(u0));
plot(p1(:,1),p1(:,2),'.-g','LineWidth', 2, 'MarkerSize', 20)
plot(p1(end,1),p1(end,2),'or','LineWidth', 2, 'MarkerSize', 15)

% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>