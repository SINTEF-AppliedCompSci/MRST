% Setting up and solving the 2D flow vertical rectangular system (2dfvrs).
% In MATLAB, this file produces Figure 8 in [A]. In GNU Octave, this file
% creates and prints the results in the folder vtk_micp_2Dfvrs which can
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
L = 500;                   % Aquifer length, m
H = 30;                    % Aquifer heigth, m
G = tensorGrid([0:.5:200 210:10:L],0:.5:H,[0 1]);
G = computeGeometry(G);
c = G.cells.centroids;
C = ones(G.cells.num,1);

% Rock
K0 = 2e-14*C;                % Aquifer permeability, m^2
porosity = .15;              % Aquifer porosity, [-]
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
fluid.crit = .1;             % Critical porosity, [-]
fluid.kmin = 1e-20;          % Minimum permeability, m^2
fluid.cells = C;             % Array with all cells, [-]
fluid.ptol = 1e-4;           % Porosity tolerance to stop the simulation

% Porosity-permeability relationship
fluid.K = @(poro) (K0.*((poro-fluid.crit)/(porosity-fluid.crit))...
        .^fluid.eta+fluid.kmin).*K0./(K0+fluid.kmin).*(poro>fluid.crit)+...
                                            fluid.kmin.*(poro<=fluid.crit);

% Gravity
gravity on
gravity y

% Strategy A (for figure a)

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f,1))>eps & G.faces.centroids(f,1)>L-.01);
fp = G.faces.centroids(f,2) * fluid.rhoWS * norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.m = zeros(size(bc.sat,1), 1);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Create Well
Q = 5e-3;     % Injection rate m^3/s
r = 0.15;     % Well radius, m
cellsW =  1:1:G.cells.num;
cellsW = cellsW(c(:,1)<1.1*c(1,1) & c(:,2)<H/10);
W = addWell([], G, rock, cellsW, 'Type', 'rate', 'Comp_i', [1,0], ...
                                           'Val', Q, 'Radius',r,'dir','y');
W.m = 0.01;
W.o = 0;
W.u = 0;
G.injectionwellonboundary = 1;
G.cellsinjectionwell = cellsW; 
                                       
% Setup some schedule
dt = hour;
nt = 400*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);
                                      
% Well different rates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are:time, rate, m, o, u.
M(1,1) = 15*hour/dt; 
M(1,2) = Q;
M(2,1) = 26*hour/dt; 
M(2,2) = eps; 
M(3,1) = 100*hour/dt; 
M(3,2) = Q;
M(3,4) = 0.04;
M(4,1) = 130*hour/dt;
M(4,2) = Q;
M(5,1) = 135*hour/dt; 
M(5,2) = eps; 
M(6,1) = 160*hour/dt; 
M(6,2) = Q;
M(6,5) = 300;
M(7,1) = 200*hour/dt; 
M(7,2) = Q;
M(8,1) = 210*hour/dt; 
M(8,2) = eps; 

% Maximum injected oxygen and urea concentrations.
fluid.Comax = max(M(:, 4));             
fluid.Cumax = max(M(:, 5));

% Create model
model = MICPModel(G, rock, fluid);

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

% Initial condition
state0 = initState(G, W, c(:,2) * fluid.rhoWS * norm(gravity), [1, 0]);
state0.m = zeros(G.cells.num,1);
state0.o = zeros(G.cells.num,1);
state0.u = zeros(G.cells.num,1);
state0.b = zeros(G.cells.num,1);
state0.c = zeros(G.cells.num,1);

% Simulate case (GNU Octave/MATLAB)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 270);
end
[~, stateA] = simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);
statea=stateA{end};

% Strategy B (for figure b)

% Create well
Whu = 1/10;     
Whb = 1 - Whu; 
cellsW =  1:1:G.cells.num;
cellsWu = cellsW(c(:,1)<1.1*c(1,1) & c(:,2) < Whu*H);
W = addWell([], G, rock, cellsWu, 'Type', 'rate', 'Comp_i', [1,0], ...
                                     'Val', Whu*Q, 'Radius', r,'dir','y');
cellsWb = cellsW(c(:,1)<1.1*c(1,1) & c(:,2) > Whu*H);
W = addWell(W, G, rock, cellsWb, 'Type', 'rate', 'Comp_i', [1,0], ...
                                     'Val', Whb*Q, 'Radius', r,'dir','y');

for i=1:2
    W(i).o = 0;
    W(i).u = 0;
    W(i).m = 0;
end
W(1).m = 0.01;                                 
G.injectionwellonboundary = 1;
G.cellsinjectionwell = [cellsWu cellsWb]; 

% Create model
model = MICPModel(G, rock, fluid);

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1) = schedule.control(i);
    schedule.control(i+1).W(1).val = Whu*M(i,2);
    schedule.control(i+1).W(2).val = Whb*M(i,2);
    schedule.control(i+1).W(1).m = M(i,3);
    schedule.control(i+1).W(1).o = M(i,4);
    schedule.control(i+1).W(1).u = M(i,5);
    schedule.step.control(M(i,1):end) = i+1;
end    

% Initial condition
state0   = initState(G, W, c(:,2) * fluid.rhoWS * norm(gravity), [1, 0]);
state0.m = zeros(G.cells.num,1);
state0.o = zeros(G.cells.num,1);
state0.u = zeros(G.cells.num,1);
state0.b = zeros(G.cells.num,1);
state0.c = zeros(G.cells.num,1);

% Simulate case (GNU Octave/MATLAB)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ok = 'true';
    fn = checkCloggingMICP(ok);
else
    fn = getPlotAfterStepMICP(state0, model, 0, 270);
end
[~, stateB] = simulateScheduleAD(state0, model, schedule,'afterStepFn',fn);
stateb=stateB{end};

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_2Dfvrs;
    cd vtk_micp_2Dfvrs;
    %mrsttovtk(G,stateA,'stateA','%f');
    mrsttovtk(G,stateB,'stateB','%f');
    return
end

% Figure 8 paper (MATLAB)
figure;
ccc=flipud(jet);
ccc=ccc(70:1:end,:);
set(gcf,'PaperUnits','inches','PaperSize',[6.83 1.85],'PaperPosition', ...
                                                          [0 0 6.83 1.85]);
set(gca,'FontName','Arial');
n1=subplot(1,2,1);
view(0, 270);
colormap (n1,ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.48 .25 .01 .5],'YTick',[0 25 50 75 100]);
xlabel({'x [m]'; '(a)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,100*(1-fluid.K(porosity-statea.c-statea.b)./K0));
s.EdgeColor = 'none';
title('Permeability reduction (using strategy A)','FontSize',8, ...
                                 'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'color','none','FontName','Arial');
zlim([0,H]);
line([90 90], [0 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([110 110], [0 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([90 110], [H/10 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([90 110], [0 0], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
n2=subplot(1,2,2);
view(0, 270);
colormap (n2,ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.93 .25 .01 .5],'YTick',[0 25 50 75 100]);
xlabel({'x [m]'; '(b)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,100*(1-fluid.K(porosity-stateb.c-stateb.b)./K0));
s.EdgeColor = 'none';
title('Permeability reduction (using strategy B)','FontSize',8, ...
                                 'FontName','Arial','Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'color','none','FontName','Arial');
zlim([0,H]);
line([90 90], [0 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([110 110], [0 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([90 110], [H/10 H/10], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
line([90 110], [0 0], [1 1],'Color','[0 0 0]','LineStyle','-', ...
                                                            'LineWidth',3);
%print -depsc2 Fig8.eps