mrstModule add co2lab ad-props;

gravity on
rhoW = 1100; 
dz = 600; 
z = linspace(300, 3000, 100)'; 
dz = 600; % max(Gt.cells.z) - min(Gt.cells.z)
Temp = z * 30 / 1e3 + 274 + 12; 
p0 = z(:) * norm(gravity) * rhoW; 
rhoGS = 760; 
bG = boCO2(Temp, rhoGS); 
figure(44), clf, 
rhoG = bG(p0) * rhoGS; 
plot(rhoG, z, 'k'); axis([0 1000 0 3000]); 
hold on; 
frhoG = @(zz) interpTable(z, rhoG, zz); 
zz = linspace(2300 - 600, 2300, 100)'; 
plot(frhoG(zz), zz, 'b', 'LineWidth', 2)
vec = [2300 - 600, 2300]; 

plot(frhoG(vec), vec, 'r*')
zz = linspace(1300 - 600, 1300, 100)'; 
plot(frhoG(zz), zz, 'r', 'LineWidth', 2)
vec = [1300 - 600, 1300]; 
plot(frhoG(vec), vec, 'b*')

set(gca, 'FontSize', 16)
set(gca, 'YDir', 'reverse')
xlabel('Density (kg / m^3)', 'FontSize', 16)
ylabel('Depth (m)','FontSize',16)
print('-depsc2','density_1D_example.eps')