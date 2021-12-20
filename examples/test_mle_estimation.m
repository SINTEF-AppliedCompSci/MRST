% Use maximum likelihood estimation to determine optimal parameters of a 
% probability distribution.  Whereas any probability distribution with any
% type and number of parameters could in principle be used, we here
% demonstrate mle estimation on a very simple example involving a univariate
% normal probability distribution.


%% Define the probability distribution function and its gradient wrt. parameters

% p here represents the vector (mu, sigma2)
normal_dist = @(x, p) 1/sqrt(2 * pi * p(2)) * exp( - (x - p(1)).^2/(2*p(2)));

% derivative with respect to mu
normal_dist_d_mu = @(x, p) (x-p(1)) .* exp(- (x-p(1)).^2 ./ (2 * p(2))) ./ ...
                    (sqrt(2*pi) .* p(2).^(3./2));

% derivative with respect to sigma2
normal_dist_d_sigma2 = @(x, p) ...
    exp( -(p(1) - x).^2./(2 * p(2))) .* ( (p(1)-x).^2 - p(2)) ./ ...
    (2 * sqrt(2*pi) .* p(2).^(5/2));

% 'gradient' with respect to parameters mu and sigma2
normal_dist_grad = @(x, p) [normal_dist_d_mu(x, p), normal_dist_d_sigma2(x, p)];


%% Generate samples from the distribution, for a set of 'true' parameters
mu = 1.5; % true mean
sigma2 = 0.45; % true variance

N = 500; % number of samples

samples = sqrt(sigma2) * randn(N, 1) + mu;

%% Estimating lambda using MLE

init_guess = [0, 0.1]; % [mu_init, sigma2_init];
pmin = [-2, 1e-4];
pmax = [ 2, 4];

[params, hist] = mle_prob_param(samples, normal_dist, normal_dist_grad, ...
                                init_guess, pmin, pmax);
   
%% Validating results

fprintf('Estimated mean:%f, sample mean: %f, actual mean: %f.\n', ...
        params(1), mean(samples), mu);
fprintf('Estimated variance:%f, sample variance: %f, actual variance: %f.\n', ...
        params(2), mean((samples - mean(samples)).^2), sigma2);