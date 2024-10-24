% Define the depth at which the transition occurs
transition_depth = 1210;
transition_width = 50; % Width of the transition zone

% Define the sigmoid function for smooth transition
sigmoid = @(x, x0, width) 1 ./ (1 + exp(-(x - x0) / width));

% Initialize s0
s0 = [0.2 0.8]; % Hack to enforce residual saturation

% Define z0g and z0l
z0g = [0.795    0.0723+0.0809    0.01905    0.020    0.0092    0.0014    0.0012         0.001];%[0.0368, 0.8430, 0.0876, 0.0243, 0.0071, 0.0006, 0.0004, 0.0002];
z0l = [0.8784    0.1023+0.0          0.0091     0.0025    0.0072    0.0004    0.0002         0];%[0.975, 0.0086, 0.0040, 0.0042, 0.0002, 0.0014, 0.0022, 0.1818];
%z0 = [0.7659    0.0000    0.1727    0.0000    0.0000    0.0000    0.0614    0.0000];
% Initialize z0 with z0l
z0 = repmat(z0l, model.G.cells.num, 1);

% Calculate the transition weights
depths = model.G.cells.centroids(:, 3);
transition_weights = 1.0- sigmoid(depths, transition_depth, transition_width);

% Ensure z0 respects z0l at the transition depth
transition_weights(depths == transition_depth) = 0;

% Apply the transition weights to z0
z0 = z0 .* (1 - transition_weights) + repmat(z0g, model.G.cells.num, 1) .* transition_weights;
% z0(:,1) = z0(:,1) + 0.35;
% z0(:,2) = z0(:,2) -0.35;
% a = z0(:,2);
% z0(:,2) = 0.5*a;
% z0(:,1) = z0(:,1) + z0(:,end) +0.5*a;
% z0(:,end) = 0.*z0(:,end) +1.0e-8;
% Display the transition weights for verification

% Now you can use z0 in your optimization or further calculations
