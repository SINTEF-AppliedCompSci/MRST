% Setting up and solving the 1D flow horizontal system (1Dfhs).
% This file produces Figure 5 in the publication.
%
% The example assumes MRST is the Matlab path. For information on
% MRST-functions, confer the MRST documentation at
%   http://www.sintef.no/projectweb/mrst/
%
%{
Copyright 2020, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling. 

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

% Required modules
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui

clear

% Setup some schedule
dt = 1*minute;
nt = 500*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Grid 
L = 75;
H = 1;
G = tensorGrid([0 .5 .55:.05:L], [0 1], [0 1]);
G = computeGeometry(G);
C = ones(G.cells.num,1);

% Rock
K0 = 1e-12*C;
porosity = 0.2;
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s                            
fluid.bW   =  @(p) 0*p + 1;  % Water formation volume factor, [-]
fluid.bO   =  @(p) 0*p + 1;  % CO2 formation volume factor, [-]
fluid.rhoWS = 1045;          % Water density, kg/m^3
fluid.rhoOS = 479;           % CO2 density, kg/m^3

% Remaining model parameters (we put them on the fluid structure)
fluid.rho_b = 35;            % Density (biofilm), kg/m^3
fluid.rho_c = 2710;          % Density (calcite), kg/m^3
fluid.k_str = 2.6e-10;       % Detachment rate, m/(Pa s)
fluid.diffm = 2.1e-9;        % Diffusion coefficient (microbes), m^2/s
fluid.diffo = 2.32e-9;       % Diffusion coefficient (oxygen), m^2/s
fluid.diffu = 1.38e-9;       % Diffusion coefficient (urea), m^2/s
fluid.alphaL = 1e-3;         % Disperison coefficient (longitudinal), m
fluid.alphaT = 4e-4;         % Disperison coefficient (transverse), m
fluid.eta = 3;               % Fitting factor, [-]
fluid.k_o = 2e-5;            % Half-velocity constant (oxygen), kg/m^3
fluid.k_u = 21.3;            % Half-velocity constant (urea), kg/m^3
fluid.mu = 4.17e-5;          % Maximum specific growth rate, 1/s
fluid.mu_u = 0.0161;         % Maximum rate of urease utilization, 1/s
fluid.k_a = 8.37e-8;         % Microbial attachment rate, 1/s                                         
fluid.k_d = 3.18e-7;         % Microbial death rate, 1/s
fluid.Y = 0.5;               % Yield growth coefficient, [-]
fluid.Yuc = 1.67;            % Yield coeccifient (calcite/urea), [-]
fluid.F = 0.5;               % Oxygen consumption factor, [-]
fluid.crit = .1;             % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation  

% Porosity-permeability relationship
fluid.K = @(poro) (K0.*((poro-fluid.crit)/(porosity-fluid.crit))...
        .^fluid.eta+fluid.kmin).*K0./(K0+fluid.kmin).*(poro>fluid.crit)+...
                                            fluid.kmin.*(poro<=fluid.crit);

% Maximum values (to ease the convergence of the solution)
fluid.omax = .04;                 % Maximum injected oxygen concentration 
fluid.umax = 300;                 % Maximum injected urea concentration
fluid.mmax = 105;                 % Maximum value of biomass concentration
fluid.bmax = porosity-fluid.ptol; % Maximum biofilm volume fraction
fluid.cmax = porosity-fluid.ptol; % Maximum calcite volume fraction 

% Create well
Q1 = 2.4e-05; % Injection rate m^3/s
Q2 = eps;     % Closed well
Cm = 0.01;    % Injected microbial concentration kg/m^3
W = addWellMICP([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1,0],...
    'Val', Q1,'o', 0,'u', 0,'m', Cm,'b', 0,'c', 0, 'Radius',...
    .15*meter,'name', 'I');

% Create model
model = MICPModel(G, rock, fluid);

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f,1))>eps & G.faces.centroids(f,1)>L-1);
fp = G.faces.centroids(f,2)*fluid.rhoWS*norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Well different rates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are: time, rate, o, u, m.
M(1,1) = 20*hour; 
M(1,2) = Q1;
M(2,1) = 40*hour; 
M(2,2) = Q2; 
M(3,1) = 140*hour; 
M(3,2) = Q1;
M(3,3) = fluid.omax;
M(4,1) = 160*hour;
M(4,2) = Q1;
M(5,1) = 180*hour; 
M(5,2) = Q2; 
M(6,1) = 230*hour; 
M(6,2) = Q1;
M(6,4) = fluid.umax;
M(7,1) = 250*hour; 
M(7,2) = Q1;
M(8,1) = 270*hour; 
M(8,2) = Q2; 

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1)=schedule.control(i);
    schedule.control(i+1).W(1).val=M(i,2);
    schedule.control(i+1).W(1).sign=sign(M(i,2));
    schedule.control(i+1).W(1).o=M(i,3);
    schedule.control(i+1).W(1).u=M(i,4);
    schedule.control(i+1).W(1).m=M(i,5);
    schedule.step.control(M(i,1)/dt:end)=i+1;
end    

% Initial condition
state0      = initState(G, W, ...
                 G.cells.centroids(:,3)*fluid.rhoWS*norm(gravity), [1, 0]);
state0.o    = zeros(G.cells.num,1);
state0.u    = zeros(G.cells.num,1);
state0.m    = zeros(G.cells.num,1);
state0.b    = zeros(G.cells.num,1);
state0.c    = zeros(G.cells.num,1);

% Simulate
[~, states] = simulateScheduleADMICP(state0, model, schedule);

% Figure 5 paper
set(gcf,'PaperUnits','inches','PaperSize',[6.83 6],'PaperPosition', ...
                                                             [0 0 6.83 6]);
set(gca,'FontName','Arial','FontSize',9);
n1=subplot(3,3,4);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)/dt}.m,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt}.m,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt}.m,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt}.m,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt}.m,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt}.m,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.m,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 .01],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 .01],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [.01 .01],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .01]);
xlabel({'x [m]';'(a)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_m$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Microbes','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                            (0:.002:.01));
n2=subplot(3,3,5);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)/dt}.o,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt}.o,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt}.o,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt}.o,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt}.o,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt}.o,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.o,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 fluid.omax],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 fluid.omax],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [fluid.omax fluid.omax],'Color','red','LineStyle',':', ...
                                                            'LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 fluid.omax]);
xlabel({'x [m]';'(b)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_o$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Oxygen','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                       (0:.01:fluid.omax));
n3=subplot(3,3,6);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)/dt}.u,'color',[0 .8 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt}.u,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt}.u,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt}.u,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt}.u,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt}.u,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt}.u,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt}.u,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.u,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 fluid.umax],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 fluid.umax],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [fluid.umax fluid.umax],'Color','red','LineStyle',':', ...
                                                            'LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 fluid.umax]);
xlabel({'x [m]';'(c)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_u$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Urea','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                        (0:60:fluid.umax));
n4=subplot(3,3,7);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)/dt}.b,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt}.b,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt}.b,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt}.b,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt}.b,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt}.b,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt}.b,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt}.b,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.b,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 .003],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 .003],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [.003 .003],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .0003]);
xlabel({'x [m]';'(d)'},'FontSize',9,'Interpreter','latex');        
ylabel('$\phi_b$ [$-$]','FontSize',9,'Interpreter','latex');
grid on
title('Biofilm','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                          (0:.0001:.0003));
n5=subplot(3,3,8);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)/dt}.c,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt}.c,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt}.c,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt}.c,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt}.c,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt}.c,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt}.c,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt}.c,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.c,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 0.04],'Color','red','LineStyle',':','LineWidth',2);                                                  
line([15 15], [0 0.04],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0.04 0.04],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .04]);
xlabel({'x [m]';'(e)'}','FontSize',9,'Interpreter','latex');        
ylabel('$\phi_c$ [$-$]','FontSize',9,'Interpreter','latex');
grid on
title('Calcite','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                              (0:.01:.05));                                                        
n6=subplot(3,3,9);
hold on
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(1,1)/dt}.b...
                            -states{M(1,1)/dt}.c)./K0),'color',[0 .8 0],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(2,1)/dt}.b...
                           -states{M(2,1)/dt}.c)./K0),'color',[0 .74 1],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(3,1)/dt}.b...
                             -states{M(3,1)/dt}.c)./K0),'color',[0 0 0],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(4,1)/dt}.b...
                           -states{M(4,1)/dt}.c)./K0),'color',[1 .5 .9],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(5,1)/dt}.b...
                           -states{M(5,1)/dt}.c)./K0),'color',[0 .74 1],...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(6,1)/dt}.b...
                            -states{M(6,1)/dt}.c)./K0),'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(7,1)/dt}.b...
                           -states{M(7,1)/dt}.c)./K0),'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{M(8,1)/dt}.b...
                          -states{M(8,1)/dt}.c)./K0),'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{end}.b...
                                   -states{end}.c)./K0),'color',[0 0 0],...
                                           'LineWidth',2,'LineStyle','--');
line([10 10], [0 100],'Color','red','LineStyle',':','LineWidth',2);                                                  
line([15 15], [0 100],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [100 100],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 100]);
xlabel({'x [m]';'(f)'}','FontSize',9,'Interpreter','latex');        
ylabel('$|\Delta K/K_0|$ [$\%$]','FontSize',9,'Interpreter','latex');
grid on
title('Permeability','FontSize',9,'FontName','Arial','Interpreter','latex')
cb=legend('$t^I_1=\;\;20\;$h','$t^I_2=\;\;40\;$h', ...
               '$t^I_3=140\;$h$\qquad \qquad \qquad$','$t^I_4=160\;$h', ...
               '$t^I_5=180\;$h','$t^I_6=230\;$h$\qquad \qquad \qquad$', ...
                    '$t^I_7=250\;$h','$t^I_8=270\;$h','$t^I_9=500\;$h', ...
                     'Location','best','Interpreter','latex','FontSize',9);
set(cb,'position',[.5 .67 .01 .15]);
cb.NumColumns = 3;
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick',...
                                                               (0:20:100));
%print -depsc2 Fig5.eps
