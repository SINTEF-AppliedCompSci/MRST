% Setting up and solving the system to compare to OPM.
% In MATLAB, this file produces Figure 3 in [B]. In GNU Octave, this file
% creates and prints the results in the folder vtk_micp_mrst_opm which can
% be visualized using ParaView.
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
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui

clear

% Setup some schedule
dt = 1*hour;
nt = 500*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Grid 
L = 100;
H = 1;
G = tensorGrid(0:L, [0 1], [0 1]);
G = computeGeometry(G);
C = ones(G.cells.num,1);

% Rock
K0 = 1e-14*C;
porosity = 0.15;
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
fluid.diffm = 0;             % Diffusion coefficient (microbes), m^2/s
fluid.diffo = 0;             % Diffusion coefficient (oxygen), m^2/s
fluid.diffu = 0;             % Diffusion coefficient (urea), m^2/s
fluid.alphaL = 0;            % Disperison coefficient (longitudinal), m
fluid.alphaT = 0;            % Disperison coefficient (transverse), m
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
Q1 = 2/day;   % Injection rate m^3/s
Q2 = eps;     % Closed well
Cm = .01;     % Injected microbial concentration kg/m^3
W = addWellMICP([], G, rock, 1, 'Type', 'rate', 'Comp_i', [1,0],...
    'Val', Q1,'o', 0,'u', 0,'m', Cm,'b', 0,'c', 0, 'Radius',...
    .15*meter,'name', 'I');
W = addWellMICP(W, G, rock, G.cells.num, 'Type', 'bhp', 'Comp_i', [1,0],...
    'Val', atm,'o', 0,'u', 0,'m', 0,'b', 0,'c', 0, 'Radius',...
    .15*meter,'name', 'P');

% Create model
model = MICPModel(G, rock, fluid);

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
schedule = simpleSchedule(timesteps,'W',W);
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
state0      = initState(G, W, atm, [1, 0]);
state0.o    = zeros(G.cells.num,1);
state0.u    = zeros(G.cells.num,1);
state0.m    = zeros(G.cells.num,1);
state0.b    = zeros(G.cells.num,1);
state0.c    = zeros(G.cells.num,1);

[~, states] = simulateScheduleADMICP(state0, model, schedule);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_mrst_opm;
    cd vtk_micp_mrst_opm;
    mrsttovtk(G,states,'states','%f');
    return
end

% Read vtk data to states_OPM from the OPM simulation (MATLAB)
fid = fopen('micp_opm_vtk/MICP_OPM.pvd', 'r');
c1 = fscanf(fid,'%c');
fclose(fid);
newStr = extractBetween(c1,'<DataSet timestep="','" file');
outNums = str2double(newStr);
for i=0:max(size(outNums))-1
    filename = sprintf('micp_opm_vtk/MICP_OPM-%05d.vtu',i);
    fileID = fopen(filename,'r');
    c1 = fscanf(fileID,'%c');
    newStr = extractBetween(c1,['<DataArray type="Float32" Name' ...
        '="urea concentration" NumberOfComponents="1" format="ascii">'],...
                                                           '</DataArray>');
    states_OPM{i+1}.u = sscanf(newStr{1},'%f');
    newStr = extractBetween(c1,['<DataArray type="Float32" Name' ...
    '="bacteria concentration" NumberOfComponents="1" format="ascii">'],...
                                                           '</DataArray>');
    states_OPM{i+1}.m = sscanf(newStr{1},'%f');
    newStr = extractBetween(c1,['<DataArray type="Float32" Name' ...
     '="oxygen concentration" NumberOfComponents="1" format="ascii">'], ...
                                                           '</DataArray>');
    states_OPM{i+1}.o = sscanf(newStr{1},'%f');
    newStr = extractBetween(c1,['<DataArray type="Float32" Name' ...
      '="biofilm" NumberOfComponents="1" format="ascii">'],'</DataArray>');
    states_OPM{i+1}.b = sscanf(newStr{1},'%f');
    newStr = extractBetween(c1,['<DataArray type="Float32" Name' ...
      '="calcite" NumberOfComponents="1" format="ascii">'],'</DataArray>');
    states_OPM{i+1}.c = sscanf(newStr{1},'%f');
    fclose(fileID);
end

% Figure (MATLAB)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9.11 12]);
set(gca,'FontName','Arial','FontSize',9);
n1=subplot(4,2,2);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),states_OPM{i}.m(1:2:end,1),...
          'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),states_OPM{end}.m(1:2:end,1),'o',...
              'MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),states{M(1,1)/dt-1}.m,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');                                        
plot(G.cells.centroids(:,1),states{M(2,1)/dt-1}.m,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt-1}.m,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt-1}.m,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt-1}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt-1}.m,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt-1}.m,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt-1}.m,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.m,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([12.5 12.5], [0 .01],'Color','red','LineStyle',':','LineWidth',2);
line([17.5 17.5], [0 .01],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [.01 .01],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .01]);
xlabel({'x [m]';'(a)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_m$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
ax = gca;
ax.YAxis.Exponent = -2;
grid on
title('Microbes','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                            (0:.002:.01));
n2=subplot(4,2,3);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),states_OPM{i}.o(1:2:end,1),...
          'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),states_OPM{end}.o(1:2:end,1),'o',...
              'MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),states{M(1,1)/dt-1}.o,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');                                        
plot(G.cells.centroids(:,1),states{M(2,1)/dt-1}.o,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');                                      
plot(G.cells.centroids(:,1),states{M(3,1)/dt-1}.o,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':'); 
plot(G.cells.centroids(:,1),states{M(4,1)/dt-1}.o,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');                                        
plot(G.cells.centroids(:,1),states{M(5,1)/dt-1}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');                                         
plot(G.cells.centroids(:,1),states{M(6,1)/dt-1}.o,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');                                       
plot(G.cells.centroids(:,1),states{M(7,1)/dt-1}.o,'color',[1 .9 0], ... 
                                            'LineWidth',2,'LineStyle','-');                                          
plot(G.cells.centroids(:,1),states{M(8,1)/dt-1}.o,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');                                        
plot(G.cells.centroids(:,1),states{end}.o,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');                                               
line([12.5 12.5], [0 fluid.omax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([17.5 17.5], [0 fluid.omax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([12.5 17.5], [fluid.omax fluid.omax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 fluid.omax]);
xlabel({'x [m]';'(b)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_o$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
ax = gca;
ax.YAxis.Exponent = -2;
grid on
title('Oxygen','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                       (0:.01:fluid.omax));
n3=subplot(4,2,4);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),states_OPM{i}.u(1:2:end,1),...
          'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),states_OPM{end}.u(1:2:end,1),'o',...
              'MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),states{M(1,1)/dt-1}.u,'color',[0 .8 0], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt-1}.u,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt-1}.u,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt-1}.u,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt-1}.u,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt-1}.u,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt-1}.u,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt-1}.u,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.u,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');                                                  
line([12.5 12.5], [0 fluid.umax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([17.5 17.5], [0 fluid.umax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([12.5 17.5], [fluid.umax fluid.umax],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 fluid.umax]);
xlabel({'x [m]';'(c)'},'FontSize',9,'Interpreter','latex');        
ylabel('$c_u$ [kg/m$^3$]','FontSize',9,'Interpreter','latex');
grid on
title('Urea','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                        (0:60:fluid.umax));
n4=subplot(4,2,5);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),states_OPM{i}.b(1:2:end,1),...
          'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),states_OPM{end}.b(1:2:end,1),'o',...
              'MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),states{M(1,1)/dt-1}.b,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt-1}.b,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt-1}.b,'color',[0 0 0], ... 
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt-1}.b,'color',[1 .5 .9], ... 
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt-1}.b,'color',[0 .74 1], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt-1}.b,'color',[0 0 0], ... 
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt-1}.b,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt-1}.b,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.b,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([12.5 12.5], [0 .00015],'Color','red','LineStyle',':','LineWidth',2);
line([17.5 17.5], [0 .00015],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [.00015 .00015],'Color','red','LineStyle',':',...
                                                            'LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .00015]);
xlabel({'x [m]';'(d)'},'FontSize',9,'Interpreter','latex');        
ylabel('$\phi_b$ [$-$]','FontSize',9,'Interpreter','latex');
grid on
title('Biofilm','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                        (0:.00003:.00015));
n5=subplot(4,2,6);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),states_OPM{i}.c(1:2:end,1),...
          'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),states_OPM{end}.c(1:2:end,1),'o',...
              'MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),states{M(1,1)/dt-1}.c,'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(2,1)/dt-1}.c,'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(3,1)/dt-1}.c,'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),states{M(4,1)/dt-1}.c,'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(5,1)/dt-1}.c,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(6,1)/dt-1}.c,'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),states{M(7,1)/dt-1}.c,'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),states{M(8,1)/dt-1}.c,'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),states{end}.c,'color',[0 0 0],'LineWidth', ...
                                                       2,'LineStyle','--');
line([12.5 12.5], [0 0.02],'Color','red','LineStyle',':','LineWidth',2);                                                  
line([17.5 17.5], [0 0.02],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [0.02 0.02],'Color','red','LineStyle',':', ...
                                                            'LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 .02]);
xlabel({'x [m]';'(e)'}','FontSize',9,'Interpreter','latex');        
ylabel('$\phi_c$ [$-$]','FontSize',9,'Interpreter','latex');
ax = gca;
ax.YAxis.Exponent = -2;
ax.YAxis.TickLabelFormat = '%.1f';
grid on
title('Calcite','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                            (0:.005:.02));                                                        
n7=subplot(4,2,7);
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-');
for i=2:9
    plot(G.cells.centroids(1:2:end,1),100/porosity*(states_OPM{i}.c ...
             (1:2:end,1)+states_OPM{i}.b(1:2:end,1)),'o','MarkerSize',5,...
                             'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
plot(G.cells.centroids(1:2:end,1),100/porosity*(states_OPM{end}.c ...
           (1:2:end,1)+states_OPM{end}.b(1:2:end,1)),'o','MarkerSize',5,...
                             'MarkerEdgeColor',[0 0 1],'LineStyle','none');
plot(G.cells.centroids(:,1),100/porosity*(states{M(1,1)/dt-1}.c+...
                               states{M(1,1)/dt-1}.b),'color',[0 .8 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100/porosity*(states{M(2,1)/dt-1}.c+...
                              states{M(2,1)/dt-1}.b),'color',[0 .74 1], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100/porosity*(states{M(3,1)/dt-1}.c+...
                                states{M(3,1)/dt-1}.b),'color',[0 0 0], ...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100/porosity*(states{M(4,1)/dt-1}.c+...
                              states{M(4,1)/dt-1}.b),'color',[1 .5 .9], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100/porosity*(states{M(5,1)/dt-1}.c+...
                              states{M(5,1)/dt-1}.b),'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100/porosity*(states{M(6,1)/dt-1}.c+...
                                states{M(6,1)/dt-1}.b),'color',[0 0 0], ...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100/porosity*(states{M(7,1)/dt-1}.c+...
                               states{M(7,1)/dt-1}.b),'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100/porosity*(states{M(8,1)/dt-1}.c+...
                              states{M(8,1)/dt-1}.b),'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),100/porosity*(states{end}.c+states{end}.b),...
                          'color',[0 0 0],'LineWidth', 2,'LineStyle','--');
line([12.5 12.5], [0 20],'Color','red','LineStyle',':','LineWidth',2);                                                  
line([17.5 17.5], [0 20],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [20 20],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 20]);
xlabel({'x [m]';'(e)'}','FontSize',9,'Interpreter','latex');        
ylabel('$\phi$ [$\%$]','FontSize',9,'Interpreter','latex');
grid on
title('Porosity','FontSize',9,'FontName','Arial','Interpreter','latex');
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick', ...
                                                              (0:5:20)); 
n8=subplot(4,2,8);
hold on                                     
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
     states{M(1,1)/dt-1}.b-states{M(1,1)/dt-1}.c)./K0),'color',[0 .8 0],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
    states{M(2,1)/dt-1}.b-states{M(2,1)/dt-1}.c)./K0),'color',[0 .74 1],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
      states{M(3,1)/dt-1}.b-states{M(3,1)/dt-1}.c)./K0),'color',[0 0 0],...
                                            'LineWidth',2,'LineStyle',':');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
    states{M(4,1)/dt-1}.b-states{M(4,1)/dt-1}.c)./K0),'color',[1 .5 .9],...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
    states{M(5,1)/dt-1}.b-states{M(5,1)/dt-1}.c)./K0),'color',[0 .74 1],...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
      states{M(6,1)/dt-1}.b-states{M(6,1)/dt-1}.c)./K0),'color',[0 0 0],...
                                           'LineWidth',2,'LineStyle','-.');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
    states{M(7,1)/dt-1}.b-states{M(7,1)/dt-1}.c)./K0),'color',[1 .9 0], ...
                                            'LineWidth',2,'LineStyle','-');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-...
   states{M(8,1)/dt-1}.b-states{M(8,1)/dt-1}.c)./K0),'color',[0 .74 1], ...
                                           'LineWidth',2,'LineStyle','--');
plot(G.cells.centroids(:,1),100*(1-fluid.K(porosity-states{end}.b...
      -states{end}.c)./K0),'color',[0 0 0],'LineWidth',2,'LineStyle','--');
for i=2:9
    k9=100*(1-fluid.K(porosity-states_OPM{i}.b-...
                                               states_OPM{i}.c)./K0);  
    plot(G.cells.centroids(1:2:end,1),k9(1:2:end),'o','MarkerSize',5,...
                             'MarkerEdgeColor',[0 0 1],'LineStyle','none');
end
k9=100*(1-fluid.K(porosity-states_OPM{end}.b-states_OPM{end}.c)./K0); 
plot(G.cells.centroids(1:2:end,1),k9(1:2:end),'o','MarkerSize',5,...
                             'MarkerEdgeColor',[0 0 1],'LineStyle','none');
line([12.5 12.5], [0 100],'Color','red','LineStyle',':','LineWidth',2);                                                  
line([17.5 17.5], [0 100],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [100 100],'Color','red','LineStyle',':','LineWidth',2);
line([12.5 17.5], [0 0],'Color','red','LineStyle',':','LineWidth',2);
hold off
xlim([0 L]);
ylim([0 100]);
xlabel({'x [m]';'(f)'}','FontSize',9,'Interpreter','latex');        
ylabel('$|\Delta K/K_0|$ [$\%$]','FontSize',9,'Interpreter','latex');
grid on
title('Permeability','FontSize',9,'FontName','Arial','Interpreter','latex')
cb=legend('$20\;$h','$40\;$h', ...
               '$140\;$h$\qquad \qquad \qquad$','$160\;$h', ...
               '$180\;$h','$230\;$h$\qquad \qquad \qquad$', ...
                    '$250\;$h','$270\;$h','$500\;$h', ...
                     'Location','best','Interpreter','latex','FontSize',9);
set(cb,'position',[.14 .78 .3 .15]);
cb.NumColumns = 3;
set(gca,'FontSize',9,'FontName','Arial','XTick',(0:20:L),'YTick',...
                                                               (0:20:100));
ax8=axes('position',[.89 .258 .01 .01],'XColor', 'none','YColor','none');
hold on
plot(0,0,'color',[0 0 0],'LineStyle','-'); 
plot(-10,-200,'o','MarkerSize',5,'MarkerEdgeColor',[0 0 1],'LineStyle',...
                                                                   'none');
hold off
legend('MRST','OPM','Location','northeast');
set(gca,'FontSize',9,'FontName','Times','XTick',(0:20:L),'YTick', ...
                                                                 (0:5:20));
%print -depsc2 Fig3.eps    