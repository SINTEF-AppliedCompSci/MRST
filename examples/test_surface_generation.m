mrstModule add static-modeling

% In this example we demonstrate how to generate a continuous bivariate
% function (which can be regarded as a 2½D function) from a set of scattered
% sampled points, using the 'generate_trend_surface' function

%% making noisy, scattered sample points using the MATLAB logo function

% In order to generate a set of sampled points from an 'unknown' surface, we
% sample from a function readily available in MATLAB, namely the membrane
% function used in the MATLAB logo.

refsurf = membrane(1, 100);
refsurfarray = refsurf(:);

% sampling 1000 random points from the surface
num_pts = 1000;

x_ixs = ceil(rand([num_pts, 1]) * size(refsurf, 1));
y_ixs = ceil(rand([num_pts, 1]) * size(refsurf, 1));
z_val = refsurfarray(sub2ind(size(refsurf), x_ixs, y_ixs));

points = [x_ixs, y_ixs, z_val];

% We add some artificial noise, to pretend these points were obtained through
% some imprecise measurement process.  
zspan = max(z_val) - min(z_val);
noise_level = 0.02; % standard deviation 2% of full z span
noise = randn([num_pts, 1]) * noise_level * zspan;

noisy_points = points;
noisy_points(:,3) = points(:,3) + noise;

% plotting exact points (red) and noisy points (blue)
clf;
plot3(points(:,1), points(:,2), points(:,3), '.r'); hold on;
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
legend('Exact sample points', 'Noisy sample points')
set(gcf, 'position', [200 200 1000 640])

%% Generating "loose" surface from scattered points, using least squares

% The generated surface should not be forced to pass through the sample
% points, since the samples are noisy and we also care about regularity.

% Surface 'stiffness' is regulated by the 'stiffness' parameter.  The
% 'num_coefs' specify how many degrees of freedom are used to describe the
% surface in each parameter direction (control points in the underlying
% tensor product spline surface)

[num_coefs, stiffness] = deal(30, 1e-8); % this is a very low stiffness
num_samples = 200;  % We will sample the surface as a regular grid.
                    % 'num_samples' here describes the number of grid points
                    % in x and y.

% generating surface using the 'generate_trend_surface' function.  We are
% here mainly interested in the sampled grid 'zgrid', but the function also
% return the abscissae of the gridpoints, as well as an object representing
% the spline surface itself, 'splinesurf'.  The latter is an instance of a
% 'BivariateSplineFunction'. 
[xcoords, ycoords, zgrid_lowstiff, splinesurf] = ...
    generate_trend_surface(noisy_points, ...
                           num_coefs, ...
                           num_samples, stiffness);


% We plot the resulting surface.  We note that while the surface is good at
% approximating the samples, it is very irregular in shape.
figure
subplot(1, 2, 1)
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid_lowstiff', 'edgealpha', 0.1);
title('Approximating surface, low stiffness')
set(gcf, 'position', [200, 200, 1400, 500])

%% Generating "stiff" surface from scattered points, using least squares


% We create a new approximating surface, this time with a significantly
% higher stiffness value.
stiffness = 1e-5;

% generating surface
[xcoords, ycoords, zgrid_histiff, splinesurf] = ...
    generate_trend_surface(noisy_points, ...
                           num_coefs, ...
                           num_samples, stiffness);


% plot surface.  We note that the surface shape is now much more pleasing,
% and closer to the real, underlying function
subplot(1, 2, 2)
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid_histiff', 'edgealpha', 0.1);
title('Approximating surface, higher stiffness')


%% Generating surface from scattered points, using multilevel B-spline algorithm

% The 'generate_trend_surface' function is based on global optimization and
% least squares.  A less costly alternative when working with high resolution
% surfaces can be to use the heuristic MBA algorithm, as described in the
% paper [10.1109/2945.620490] (Wolberg and Shin, 1997).  We here demonstrate
% its use on the same set of noisy samples used above

box = [min(noisy_points(:,1)), min(noisy_points(:,2)), ...
       max(noisy_points(:,1)), max(noisy_points(:,2))];

mbaspline = scattered_point_approximation(noisy_points, 5, [2,2], box, false);

% The function only returns the spline object, so in order to create a grid,
% we have to sample it from the spline function ourseles
x = linspace(box(1), box(3), num_samples);
y = linspace(box(2), box(4), num_samples);
   
[U, V] = ndgrid(x, y);

zgrid = reshape(mbaspline.evaluate([U(:), V(:)], 0, 0), num_samples, num_samples);
 
% plot surface
figure
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid', 'edgealpha', 0.1);
title('Approximating surface, using MBA algorithm')

%% Copyright Notice
%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
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
