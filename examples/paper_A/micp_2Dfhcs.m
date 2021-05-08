% Setting up and solving the 2D flow horizontal circular system (2Dfhcs).
% In MATLAB, this file produces Figure 6 in [A]. In GNU Octave, this file
% creates and prints the results in the folder vtk_micp_2Dfhcs which can
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

% To get distmesh for first time, uncomment and run the following lines
% pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
% mkdir(pth)
% unzip('http://persson.berkeley.edu/distmesh/distmesh.zip', pth);
% mrstPath('reregister','distmesh', pth);

% Required modules
pth = fullfile(ROOTDIR,'utils','3rdparty','distmesh');
mrstPath('reregister','distmesh', pth);
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui ...
                                                                   distmesh
clear
                                                               
% Setup some schedule
dt = 1*hour;
nt = 500*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Grid 
L=75;
H=1;
hmin=.75;
hmid=.75;
hmax=10;
fd=@(p) dcircle(p,0,0,L);
fh=@(p) min(hmin+.3*abs(dcircle(p,0,0,0)),hmid).* ...
  (abs(dcircle(p,0,0,0))<25)+min(hmid+.3*abs(dcircle(p,0,0,25)),hmax).* ...
                                               (abs(dcircle(p,0,0,0))>=25);
[p,t]=distmesh2d(fd, fh, hmin, [-L,-L;L,L], [-L,-L;L,-L;-L,L;L,L;0,0]);
G = makeLayeredGrid(pebi(triangleGrid(p, t)),1);
close
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

% Create Well
Q1 = 1.2e-3; % Injection rate m^3/s
Q2 = eps;    % Closed well
Cm = 0.01;   % Injected microbial concentration kg/m^3
cellsWell1 =  G.cells.indexMap;
cellsWell1 = cellsWell1(G.cells.centroids(:,1).^2+...
                                       G.cells.centroids(:,2).^2<hmin^2/2);
W = addWellMICP([], G, rock, cellsWell1(1), 'Type', 'rate', 'Comp_i', ...
         [1,0], 'Val', Q1,'o', 0,'u', 0,'m', Cm,'b', 0,'c', 0, 'Radius',...
                                                    .15*meter,'name', 'I');

% Create model
model = MICPModel(G, rock, fluid);

% Boundary condition
f = boundaryFaces(G);
f = f((abs(G.faces.normals(f,2))>eps | abs(G.faces.normals(f,1))>eps) ...
                & G.faces.centroids(f,1).^2+G.faces.centroids(f,2).^2>L-1);
fp = 0;
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Well different tates and times
N = 8; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are:time, rate, o, u, m.
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
state0      = initState(G, W, G.cells.centroids(:,3)* ...
                                        fluid.rhoWS*norm(gravity), [1, 0]);
state0.o    = zeros(G.cells.num,1);
state0.u    = zeros(G.cells.num,1);
state0.m    = zeros(G.cells.num,1);
state0.b    = zeros(G.cells.num,1);
state0.c    = zeros(G.cells.num,1);

% Simulate
[~, states] = simulateScheduleADMICP(state0, model, schedule);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_2Dfhcs;
    cd vtk_micp_2Dfhcs;
    mrsttovtk(G,states,'states','%f');
    return
end

% Figure 6 paper (MATLAB)
cc(:,1)=[.75:.01:1]';
cc(:,2)=[.75:.01:1]';
cc(:,3)=[.75:-.03:0]';
ccc=flipud(jet);
ccc=ccc(70:1:end,:);
set(gcf,'PaperUnits','inches','PaperSize',[6.83 3.85],'PaperPosition', ...
                                                             [0 0 6.83 4]);
set(gca,'FontName','Arial');
n1=subplot(1,2,1);
colormap (n1,cc);
caxis([0 1e-12]);
cb = colorbar; 
title(cb, '$m^2$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.475 .25 .01 .55],'Ticks',[0 1e-12],'FontSize',8, ...
    'TickLabels',{'0','$10^{-12}$'},'TickLabelInterpreter','latex');
plotCellData(G,K0);
title('Initial permeability','Interpreter','latex','FontSize',8);
axis equal tight;
xlabel({'x [m]','{(a)}'});
ylabel('y [m]');
xlim([-75,75]);
ylim([-75,75]);
view(0,90)
set(gca,'FontSize',8,'XTick',(-75:25:75),'YTick',(-75:25:75),'color', ...
                                                'none','FontName','Arial');
n2=subplot(1,2,2);
s=plotCellData(G,100*(K0-fluid.K(porosity-states{end}.c- ...
                                                      states{end}.b))./K0);
title('Permeability (phase I MICP)','Interpreter','latex','FontSize',8);
axis equal tight;
colormap(n2,ccc)
caxis([0 100]);
cb = colorbar; 
title(cb, '$\%$','Interpreter','latex');
set(cb,'position',[.925 .25 .01 .55],'YTick',0:20:100,'FontSize',8);
xlabel({'x [m]','(b)'});
ylabel('y [m]');
s.EdgeColor = 'none';
xlim([-75,75]);
ylim([-75,75]);
view(0,90)
grid off;
set(gca,'FontSize',8,'XTick',(-75:25:75),'YTick',(-75:25:75),'color', ...
                                                'none','FontName','Arial');
rectangle('Position',[10,-2.5,5,5],'LineWidth',2,'LineStyle','-', ...
                                                   'edgecolor','[0 0 0]');
ax3=axes('position',[.15 .78 .05 .06]);
box on 
axis([-1 1 -1 1]);
xlim([-1 1])
ylim([-1 1])
colormap(ax3,'parula')
s=plotGrid(G);
view(0,90)
set(gca,'FontSize',8,'XTick',(-1:1:1),'YTick',(-1:1:1),'color', ...
                                                'none','FontName','Arial');
%print -depsc2 Fig6.eps