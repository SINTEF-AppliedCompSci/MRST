%% making noisy, scattered sample points using the MATLAB logo function

refsurf = membrane(1, 100);
refsurfarray = refsurf(:);

% sampling 1000 random points from the surface
num_pts = 1000;

x_ixs = ceil(rand([num_pts, 1]) * size(refsurf, 1));
y_ixs = ceil(rand([num_pts, 1]) * size(refsurf, 1));
z_val = refsurfarray(sub2ind(size(refsurf), x_ixs, y_ixs));

points = [x_ixs, y_ixs, z_val];

% adding noise
zspan = max(z_val) - min(z_val);
noise_level = 0.02; % standard deviation 2% of full z span
noise = randn([num_pts, 1]) * noise_level * zspan;

noisy_points = points;
noisy_points(:,3) = points(:,3) + noise;

% plotting exact points (red) and noisy points (blue)
clf;
plot3(points(:,1), points(:,2), points(:,3), '.r'); hold on;
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');

%% Generating "loose" surface from scattered points, using least squares

[num_coefs, stiffness] = deal(30, 1e-8);
num_samples = 200;

% generating surface
[xcoords, ycoords, zgrid, splinesurf] = ...
    generate_trend_surface(noisy_points, ...
                           num_coefs, ...
                           num_samples, stiffness);


% plot surface
clf;
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid', 'edgealpha', 0.1);

%% Generating "stiff" surface from scattered points, using least squares
stiffness = 1e-5;

% generating surface
[xcoords, ycoords, zgrid, splinesurf] = ...
    generate_trend_surface(noisy_points, ...
                           num_coefs, ...
                           num_samples, stiffness);


% plot surface
clf;
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid', 'edgealpha', 0.1);

%% Generating surface from scattered points, using multilevel B-spline algorithm

box = [min(noisy_points(:,1)), min(noisy_points(:,2)), ...
       max(noisy_points(:,1)), max(noisy_points(:,2))];

mbaspline = scattered_point_approximation(noisy_points, 5, [2,2], box, false);

% sampling from spline
x = linspace(box(1), box(3), num_samples);
y = linspace(box(2), box(4), num_samples);
   
[U, V] = ndgrid(x, y);

zgrid = reshape(mbaspline.evaluate([U(:), V(:)], 0, 0), num_samples, num_samples);
 
% plot surface
clf;
plot3(noisy_points(:,1), noisy_points(:,2), noisy_points(:,3), 'b.');
hold on;
surf(zgrid', 'edgealpha', 0.1);

