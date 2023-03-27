function plotRelperm(g_top,rock)
%% Generate relative permeability data.
h  = linspace(0, max(g_top.cells.H),100).';
kr = zeros([g_top.cells.num, numel(h)]);
vporo = zeros([g_top.cells.num, numel(h)]);
tic;
for m = 1 : numel(h),
   kr(:,m) = integrateVertically(rock.perm(:,1), h(m), g_top);
end
for m = 1 : numel(h),
   vporo(:,m) = integrateVertically(rock.poro, h(m), g_top);
end
toc;

%% Compute selected mobilities.
n     = 300;
krco2 = kr(1 : n : end, :);
K_av  = krco2(:,end) ./ h(end);
krw   = bsxfun(@minus, krco2(:, end), krco2);

mu_co2     = 0.1;
mu_w       = 1;
lambda_co2 = krco2 ./ mu_co2;
lambda_w   = krw   ./ mu_w;
harm       = (lambda_co2 .* lambda_w) ./ (lambda_co2 + lambda_w);

%% Plot results.
figure(1)
subplot(3,1,1)
plot(h, lambda_co2, h, lambda_w);
subplot(3,1,2)
plot(h, bsxfun(@rdivide, harm, K_av));
subplot(3,1,3)
plot(h, lambda_co2./(lambda_w+lambda_co2));

figure(2)
subplot(3,1,1)
plot(h, harm);
subplot(3,1,2)
plot(h(1:end-1),diff(vporo(1:n:end,:),1,2));
subplot(3,1,3)
plot(h(1:end-1),diff(kr(1:n:end,:),1,2));

figure(3)
plot(h,lambda_co2.*lambda_w./(lambda_w+lambda_co2));
