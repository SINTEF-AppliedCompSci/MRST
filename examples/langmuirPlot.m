clc;clear;close all;
% parameters.p_langmuir=1562*psia; %Langmuir pressure
parameters.p_langmuir=1.077e+7; % Pressure in Pa
% parameters.v_langmuir=56; %[scf/ton] Langmuir volume
parameters.v_langmuir=3; %[kg/m3] Langmuir volume

PL=parameters.p_langmuir;% Pa
VL=parameters.v_langmuir; %m3/kg

m_ad_func = @(p) (p.*VL)./(PL+p);

% p= linspace(14.5*psia, 5000*psia, 50);
p= linspace(0, 10e+7, 100);

figure(1);
size=[4,4];
% plot(convertTo(p, psia), m_ad_func(p),'bo-', 'LineWidth', 2);
plot(p, m_ad_func(p),'bo-', 'LineWidth', 2);
set(gca,'FontSize',20);
%plot(p, m_ad_func(p),'bo-', 'LineWidth', 2);plot(convertTo(p, psia), m_ad_func(p),'bo-', 'LineWidth', 2);
xlabel('Pressure, Pa')
ylabel('\rho_s, kg/m^3')
