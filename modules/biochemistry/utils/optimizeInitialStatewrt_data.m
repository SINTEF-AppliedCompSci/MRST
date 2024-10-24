% Initial guess for z0
z0_initial = [0.6745, 0.2023, 0.0905, 0.0250, 0.0072, 0.0004, 0.0002, 0]; % Replace with your initial guess


% Define the depth at which the transition occurs
transition_depth = 1210;
transition_width = 50; % Width of the transition zone

% Define the sigmoid function for smooth transition
sigmoid = @(x, x0, width) 1 ./ (1 + exp(-(x - x0) / width));

% Initialize s0
s0 = [0.2 0.8]; % Hack to enforce residual saturation

% Define z0g and z0l
z0g = [0.0368, 0.8430, 0.0876, 0.0243, 0.0071, 0.0006, 0.0004, 0.0002];
z0l = [0.7975, 0.0086, 0.0040, 0.0042, 0.0002, 0.0014, 0.0022, 0.1818];

% Initialize z0 with z0l
z0 = repmat(z0l, model.G.cells.num, 1);

% % Calculate the transition weights
% depths = model.G.cells.centroids(:, 3);
% transition_weights = 1.0- sigmoid(depths, transition_depth, transition_width);
% 
% % Ensure z0 respects z0l at the transition depth
% transition_weights(depths == transition_depth) = 0;
% 
% % Apply the transition weights to z0
% z0_initial = z0 .* (1 - transition_weights) + repmat(z0g, model.G.cells.num, 1) .* transition_weights;

% Define any constraints on z0 if necessary
% For example, if z0 should be between 0 and 1:
lb = zeros(size(z0_initial)) + 1.0e-8; % Lower bound
ub = ones(size(z0_initial));  % Upper bound

% Constraint that z0 must sum to 1
Aeq = ones(1, length(z0_initial));
beq = 1;

% Options for the optimizer
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% Run the optimization
z0_optimized = fmincon(@(z0) objectiveFunction(model, P0, T0, s0, nbact0, eos, z0, s_given), z0_initial, [], [], Aeq, beq, lb, ub, [], options);
s_given = 0.3;
% Display the optimized z0
disp('Optimized z0:');
disp(z0_optimized);

% Verify the result
state_optimized = initCompositionalStateBacteria(model, P0 * ones(model.G.cells.num, 1), T0, s0, z0_optimized, nbact0, eos);
disp('Difference between state.y and y_given:');
disp(norm(state_optimized.y - y_given));

% Define the objective function
function obj = objectiveFunction(model, P0, T0, s0, nbact0, eos, z0, s_given)
    % Initialize the state with the current z0
    state0 = initCompositionalStateBacteria(model, P0 * ones(model.G.cells.num, 1), T0, s0, z0, nbact0, eos);

    % Calculate the difference between the current state.y and y_given
    obj = norm(state0.s(:,1) - s_given);
end
