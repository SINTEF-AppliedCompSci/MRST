%% create artificial sample data
num_measurements = 3000; % number of samples
D = 6000 * meter; % distance 
sigma = 10 * meter; % correllation distance (Gaussian correllation function)                
xvals = linspace(0, D, num_measurements);

%% defining correllation function and simulating sample data

corr_fun = @(x) exp(-(x/(sigma/D)).^2); % normalize sigma by distance 
measurements = GaussianProcess1D(num_measurements, corr_fun);

%% Estimating semi-variogram
curvepts = 500;
max_reldist = 15 * sigma / D;
rel_radius = 3;
[curve, sill, range] = estimateSemivariogram1D({xvals}, {measurements}, curvepts, ...
                                           max_reldist, rel_radius);

%% Plotting results

figure;
subplot(2,1,1); 
plot(xvals, measurements); 
xlabel('x-position'); ylabel('measurement value'); title("measurements");
subplot(2,1,2);
plot(linspace(0, D * max_reldist, curvepts+1), 2 * curve); 
xlabel('distance'); ylabel('mean squared difference'); title('estimated variogram');

                                           
