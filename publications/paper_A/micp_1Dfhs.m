% Setting up and solving the 1D flow horizontal system (1Dfhs).
% In MATLAB, this file produces Figure 5 in [A]. In GNU Octave, this file
% creates and prints the results in the folder vtk_micp_1Dfhs which can
% be visualized using ParaView.
%
% The example assumes MRST is the Matlab/Octave path. For information on
% MRST-functions, confer the MRST documentation at
%   http://www.sintef.no/projectweb/mrst/
%
%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational 
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
mrstModule add ad-blackoil ad-core ad-micp

% Grid 
L = 75;                                       % Aquifer length, m
X = [0 0.5 0.55 : 0.05 : L];                  % X discretization
G = tensorGrid(X, [0 1], [0 1]);
G = computeGeometry(G);
C = ones(G.cells.num,1);

% Rock
K0 = 1e-12*C;                % Aquifer permeability, m^2
porosity = 0.2;              % Aquifer porosity, [-]
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
fluid.k_a = 8.51e-7;         % Microbial attachment rate, 1/s                                         
fluid.k_d = 3.18e-7;         % Microbial death rate, 1/s
fluid.Y = 0.5;               % Yield growth coefficient, [-]
fluid.Yuc = 1.67;            % Yield coeccifient (calcite/urea), [-]
fluid.F = 0.5;               % Oxygen consumption factor, [-]
fluid.crit = 0.1;            % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation

% Porosity-permeability relationship
fluid.K = @(poro) (K0.*((poro-fluid.crit)/(porosity-fluid.crit))...
        .^fluid.eta+fluid.kmin).*K0./(K0+fluid.kmin).*(poro>fluid.crit)+...
                                            fluid.kmin.*(poro<=fluid.crit);

% Create well
Q = 2.4e-05;  % Injection rate, m^3/s
r = 0.15;     % Well radius, m
W = addWell([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1,0], 'Val', Q, ...
                                                              'Radius', r);
W.o = 0;
W.u = 0;
W.m = 0.01;
G.injectionwellonboundary = 1;
G.cellsinjectionwell = 1; 

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f,1))>eps & G.faces.centroids(f,1) > X(end-1));
bc = addBC([], f, 'pressure', atm, 'sat', [0 0]);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Setup some schedule
dt = minute;
nt = 500*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Well different rates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are: time, rate, m, o, u.
M(1,1) = 20*hour/dt; 
M(1,2) = Q;
M(2,1) = 40*hour/dt; 
M(2,2) = eps; 
M(3,1) = 140*hour/dt; 
M(3,2) = Q;
M(3,4) = 0.04;
M(4,1) = 160*hour/dt;
M(4,2) = Q;
M(5,1) = 180*hour/dt; 
M(5,2) = eps; 
M(6,1) = 230*hour/dt; 
M(6,2) = Q;
M(6,5) = 300;
M(7,1) = 250*hour/dt; 
M(7,2) = Q;
M(8,1) = 270*hour/dt; 
M(8,2) = eps; 

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1) = schedule.control(i);
    schedule.control(i+1).W.val = M(i,2);
    schedule.control(i+1).W.m = M(i,3);
    schedule.control(i+1).W.o = M(i,4);
    schedule.control(i+1).W.u = M(i,5);
    schedule.step.control(M(i,1):end) = i+1;
end    

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 4));             
fluid.Cumax = max(M(:, 5));

% Create model
model = MICPModel(G, rock, fluid);

% Initial condition
state0      = initState(G, W, atm, [1, 0]);
state0.m    = zeros(G.cells.num,1);
state0.o    = zeros(G.cells.num,1);
state0.u    = zeros(G.cells.num,1);
state0.b    = zeros(G.cells.num,1);
state0.c    = zeros(G.cells.num,1);

% Simulate case (GNU Octave/MATLAB)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 270);
end
[~, states] = simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_1Dfhs;
    cd vtk_micp_1Dfhs;
    mrsttovtk(G,states,'states','%f');
    return
end

% Figure 5 paper (MATLAB)
figure;
set(gcf,'PaperUnits','inches','PaperSize',[6.83 6],'PaperPosition', ...
                                                             [0 0 6.83 6]);
set(gca,'FontName','Arial','FontSize',9);
n1=subplot(3,3,4);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)-1}.m,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)-1}.m,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)-1}.m,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)-1}.m,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)-1}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)-1}.m,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)-1}.m,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)-1}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.m,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 0.01],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 0.01],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0.01 0.01],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 0.01]);
xlabel({'x [m]';'(a)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_m$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Microbes','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                        (0:.002:0.01));
n2=subplot(3,3,5);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)-1}.o,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)-1}.o,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)-1}.o,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)-1}.o,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)-1}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)-1}.o,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)-1}.o,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)-1}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.o,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 0.04],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 0.04],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [0.04 0.04],'Color','red','LineStyle',':', ...
                                                            'LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 0.04]);
xlabel({'x [m]';'(b)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_o$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Oxygen','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                         (0:.01:0.04));
n3=subplot(3,3,6);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)-1}.u,'color',[0 .8 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)-1}.u,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)-1}.u,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)-1}.u,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)-1}.u,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)-1}.u,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)-1}.u,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)-1}.u,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.u,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([10 10], [0 300],'Color','red','LineStyle',':','LineWidth',2);
line([15 15], [0 300],'Color','red','LineStyle',':','LineWidth',2);
line([10 15], [300 300],'Color','red','LineStyle',':', ...
                                                            'LineWidth',2);
line([10 15], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 300]);
xlabel({'x [m]';'(c)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_u$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Urea','FontSize',9,'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:10:L),'YTick', ...
                                                          (0:60:300));
n4=subplot(3,3,7);
hold on
plot(G.cells.centroids(:,1),states{M(1,1)-1}.b,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)-1}.b,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)-1}.b,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)-1}.b,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)-1}.b,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)-1}.b,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)-1}.b,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)-1}.b,'color',[0 .74 1], ...
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
plot(G.cells.centroids(:,1),states{M(1,1)-1}.c,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)-1}.c,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)-1}.c,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)-1}.c,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)-1}.c,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)-1}.c,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)-1}.c,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)-1}.c,'color',[0 .74 1], ...
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
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(1,1)-1}.b...
                      -states{M(1,1)-1}.c)./K0)*100,'color',[0 .8 0],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(2,1)-1}.b...
                     -states{M(2,1)-1}.c)./K0)*100,'color',[0 .74 1],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(3,1)-1}.b...
                       -states{M(3,1)-1}.c)./K0)*100,'color',[0 0 0],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(4,1)-1}.b...
                     -states{M(4,1)-1}.c)./K0)*100,'color',[1 .5 .9],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(5,1)-1}.b...
                     -states{M(5,1)-1}.c)./K0)*100,'color',[0 .74 1],...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(6,1)-1}.b...
                      -states{M(6,1)-1}.c)./K0)*100,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(7,1)-1}.b...
                     -states{M(7,1)-1}.c)./K0)*100,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{M(8,1)-1}.b...
                    -states{M(8,1)-1}.c)./K0)*100,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),(1-fluid.K(porosity-states{end}.b...
                               -states{end}.c)./K0)*100,'color',[0 0 0],...
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