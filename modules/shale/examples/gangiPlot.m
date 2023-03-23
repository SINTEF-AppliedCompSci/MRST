%{
Set up matrix perm reduction due to micro-fracture clousure

Reference
---------
Gangi,1978 10.1016/0148-9062(78)90957-9
Shale data1 SPE-18267-PA
Shale data2 10.1016/j.coal.2015.12.014 table 4

Arguments
---------
Pc    -- Overburden confining pressure, Pa
m     -- Gangi surface roughness coefficient, dimensionless
P1    -- matrix effective stress when micro-fracture fully closed, Pa
alpha -- biot conefficient
p     -- pore pressure,Pa

Return
---------
[k_gangi] permeability corrfection term due to fracture clousure

Author: Bin Wang(binwang.0213@gmail.com)
Date: Dec.2018
%}


m=0.5; 
Pc=38e6; %1.1 psi/feet GEOMECHANICAL STUDIES OF THE BARNETT SHALE, TEXAS, USA
P1=180e6;
alpha=0.5;

k_gangi_func=@(p) (1-((Pc-alpha.*p)./P1).^m).^3;

p    = linspace(0, 300e5, 50);
% figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
plot(p, k_gangi_func(p),'bo-', 'LineWidth', 2)
set(gca,'FontSize',16);
grid on;
ylim([0.15 0.3])
xlabel('Pressure,Pa')
ylabel('F_{app}')
