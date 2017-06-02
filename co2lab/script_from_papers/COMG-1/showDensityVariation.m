%% Density variation with depth
% We use the PVT model of Span and Wagner to generate a plot of density
% variations with depth of the aquifer. To this end, we assume a thermal
% gradient of 30 K/km and a surface temperature of 12 centigrade. Given a
% hydrostatic pressure computed from a constant brine density of 1100
% kg/m^3, there will be different regimes of the density variation
% depending upon the depth of the aquifer
mrstModule add co2lab;

gravity on
rhoW  = 1100 * kilogram / meter^3;
z     = linspace(300, 3000, 100)' * meter;
T     = z * 30 / 1e3 + 274 + 12;
p0    = z(:) * norm(gravity) * rhoW;
co2   = SampledProp2D('rho','CarbonDioxide_100000_400000000_278_524_800_800_D');
rhoG  = co2.rho(p0, T);
frhoG = @(zz) interpTable(z, rhoG, zz);

% Plot the density variation and two intervals
z1 = 2300;
z2 = 1300;
dz = 600;
figure
plot(rhoG, z, 'k'); axis([0 1000 0 3000]);
hold on;

zz = linspace(z1 - dz, z1, 100)';
plot(frhoG(zz), zz, 'b', 'LineWidth', 2)
plot(frhoG([z1 - dz, z1]), [z1 - dz, z1], 's', 'MarkerEdgeColor', 'b',...
     'MarkerFaceColor', 'r', 'MarkerSize', 6)

zz = linspace(z2 - dz, z2, 100)';
plot(frhoG(zz), zz, 'r', 'LineWidth', 2)
plot(frhoG([z2 - dz, z2]), [z2 - dz, z2], 's', 'MarkerEdgeColor', 'r',...
     'MarkerFaceColor', 'b', 'MarkerSize', 6)
set(gca, 'YDir', 'reverse', 'FontSize', 12)
xlabel('Density [kg / m^3]', 'FontSize', 12)
ylabel('Depth [m]','FontSize',12)

% print -depsc2 density_1D_example.eps;