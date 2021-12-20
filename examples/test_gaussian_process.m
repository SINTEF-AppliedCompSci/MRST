%% Generate 1D gaussian random field
num_samples = 1000;
corr_fun = @(x) exp(-abs(x)/0.3); % we use exponential correllation function

field_1D_A = GaussianProcess1D(num_samples, corr_fun);

corr_fun = @(x) exp(-(x/0.05).^2);

field_1D_B = GaussianProcess1D(num_samples, corr_fun);
clf;
subplot(1,2,1); plot(field_1D_A); title('Exponential correllation');
subplot(1,2,2); plot(field_1D_B); title('Gaussian correllation');

%% Generate 3D gaussian random field
num_samples = 400;

% generating field
corr_fun = @(xy) exp(-sum(xy.^2, 2)/0.0001);
fieldA = GaussianProcessND([num_samples, num_samples], corr_fun);

% Re-generating field, changing correllation length
corr_fun = @(xy) exp(-sum(xy.^2, 2)/0.01);
fieldB = GaussianProcessND([num_samples, num_samples], corr_fun);

% Re-generating field using exponential correllation function
corr_fun = @(xy) exp(-sqrt(sum(xy.^2, 2))/0.1);
fieldC = GaussianProcessND([num_samples, num_samples], corr_fun);


figure;
subplot(1,3,1); surf(fieldA, 'edgealpha', 0); view(0, 90); 
title("Gaussian correllation, sigma=1e-2");
                                                  
subplot(1,3,2); surf(fieldB, 'edgealpha', 0); view(0, 90);
title("Gaussian correllation, sigma=1e-1");

subplot(1,3,3); surf(fieldC, 'edgealpha', 0); view(0, 90);
title("Exponential correllation");

