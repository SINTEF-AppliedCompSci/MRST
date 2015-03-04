% script to generate left pannel of Figure 4 of paper
% It plots the density of CO2 as function of depth with a given termal
% gradiant of 30 K/km with 12 C at the h=0 and assuming hydrostatic
% conditions
gravity on
rhoW=1100; % brine density
z=linspace(300,3000,100)'; % valuse to plot
dz=600; %max(Gt.cells.z)-min(Gt.cells.z) max change in aquifer
Temp=z*30/1e3+274+12; % temprature profile
p0=z(:)*norm(gravity)*rhoW; % hydrostatic pressure profile
%pp0 = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
rhoGS=760; % reference CO2 density used for black oil model
bG  =  boCO2(Temp, rhoGS); 

% prepear for ploting
figure(44),clf,
rhoG=bG(p0)*rhoGS;% density of CO2 in black oil formulation
plot(rhoG,z,'k');axis([0 1000 0 3000]); % plot for all values
hold on;
frhoG=@(zz) interpTable(z,rhoG,zz);
zz=linspace(2300-600,2300,100)';
% plot for values of aquifer with depth 2300
plot(frhoG(zz),zz,'b','LineWidth',2);
vec=[2300-600,2300];

plot(frhoG(vec),vec,'r*')
zz=linspace(1300-600,1300,100)';

% plot for values of aquifer with depth 2300
plot(frhoG(zz),zz,'r','LineWidth',2)
vec=[1300-600,1300];
plot(frhoG(vec),vec,'b*')

% and legends
set(gca,'FontSize',16)
set(gca,'YDir','reverse')
xlabel('Density (kg/m^3)','FontSize',16)
ylabel('Depth (m)','FontSize',16)

print('-depsc2','figs/density_1D_example.eps')