% Setting up and solving the 2D flow leaky system (2Dfls).
% In MATLAB, this file produces Figure 10 and 11 in [A]. In GNU Octave, 
% this file creates and prints the results in the folder vtk_micp_2Dfls 
% which can be visualized using ParaView.
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
mrstModule add deckformat ad-core ad-blackoil ad-micp ad-props mrst-gui

clear

% Setup some schedule
dt = 1*hour;
nt = 950*hour/dt;
clear schedule
timesteps = repmat(dt, nt, 1);

% Grid
L=500;
H=160;
if exist('OCTAVE_VERSION', 'builtin') ~= 0 %GNU Octave
    G = tensorGrid([0 1:.5:99.5 99.85 100.15 100.5:.5:203 ...
                                          L*exp(-.85:.05:0)],0:.5:H,[0 1]);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<100-.3 | c(:,1)>100+.3) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<99.9 | c(:,1)>100.1) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<-.5-eps & c(:,2)>130);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<99.8 & c(:,2)<30);
    G = computeGeometry(G);                                                       
else %MATLAB
    [X1,Y1] = meshgrid([-L:10:-50 180:10:L], [0:5:30 31:1:129 130:5:H]);
    [X2,Y2] = meshgrid(-50:10:50, 0:5:30);
    [x,y] = meshgrid([100-.3 100 100+.3], 0:.25:H);
    [xc,yc] = meshgrid([-1 0 1], 130:.25:H);
    [xw,yw] = meshgrid([-50*exp(0:-0.25:-3.6) 50*exp(-3.6:0.25:0)], ...
                                                                  130:1:H);
    [xl,yl] = meshgrid([100-50*exp(0:-0.25:-5) 100+80*exp(-5:0.125:0)], ...
                                                                    0:1:H);
    [xwc1,ywc1] = meshgrid(50*exp(-3.6:0.25:0), 130:.25:137.5);
    [xwc2,ywc2] = meshgrid(50*exp(-3.6:0.25:0), 137.5:.5:145);
    [xwc3,ywc3] = meshgrid(100-50*exp(0:-0.25:-5), 130:.25:137.5);
    [xwc4,ywc4] = meshgrid(100-50*exp(0:-0.25:-5), 137.5:.5:145);
    P = unique([X1(:) Y1(:); X2(:) Y2(:); x(:) y(:); xw(:) yw(:); xl(:) ...
          yl(:); xc(:) yc(:); xwc1(:) ywc1(:); xwc2(:) ywc2(:); xwc3(:) ...
                                        ywc3(:); xwc4(:) ywc4(:)], 'rows');
    G = triangleGrid(P);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<100-.3 | c(:,1)>100+.3) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = makeLayeredGrid(pebi(G),1);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, (c(:,1)<99.9 | c(:,1)>100.1) & (c(:,2)<130 & ...
                                                               c(:,2)>30));
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<-.5-eps & c(:,2)>130);
    G = computeGeometry(G);
    c = G.cells.centroids;
    G = removeCells(G, c(:,1)<99.8 & c(:,2)<30);
    G = computeGeometry(G);
end
C = ones(G.cells.num,1);

% Rock
K0 = 2e-14*C;
cellsfrac =  G.cells.indexMap;
cellsfrac1 = cellsfrac(G.cells.centroids(:,1)>99.9 & ... 
            G.cells.centroids(:,1)<100.1 & G.cells.centroids(:,2)<130 & ...
                                                G.cells.centroids(:,2)>30);
cellsF =  G.cells.indexMap;
idx = ismember(cellsF,cellsfrac1);
K0(idx) = 1e-12;
porosity = 0.15;
rock = makeRock(G, K0, porosity);

% Fluid properties
fluid.muw = 2.535e-4;        % Water viscocity, Pa s   
fluid.muO   = 3.95e-5;       % CO2 viscosity, Pa s
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
fluid.omax = 0.04;                % Maximum injected oxygen concentration 
fluid.umax = 300;                 % Maximum injected urea concentration
fluid.mmax = 105;                 % Maximum value of biomass concentration
fluid.bmax = porosity-fluid.ptol; % Maximum biofilm volume fraction
fluid.cmax = porosity-fluid.ptol; % Maximum calcite volume fraction

% Create well
Q1 = 6e-3;   % Injection rate m^3/s
Q2 = eps;    % Closed well
Cm = 0.01;   % Injected microbial concentration kg/m^3
neww=1:1:G.cells.num;
neww = neww(G.cells.centroids(:,1)<min(G.cells.centroids(:,1))+.1 & ...
                  G.cells.centroids(:,2)>130 & G.cells.centroids(:,2)<133);
W = addWellMICP([], G, rock, neww, 'Type', 'rate', 'Comp_i', [1,0],...
             'Val', Q1/10,'o', 0,'u', 0,'m', Cm,'b', 0,'c', 0, 'Radius',...
                                          .15*meter,'name', 'I','dir','y');
new2=1:1:G.cells.num;
new2 = new2(G.cells.centroids(:,1)<min(G.cells.centroids(:,1))+.1 & ...
                                               G.cells.centroids(:,2)>133);
W = addWellMICP(W, G, rock, new2, 'Type', 'rate', 'Comp_i', [1,0],...
            'Val', 9*Q1/10,'o', 0,'u', 0,'m', 0,'b', 0,'c', 0, 'Radius',...
                                          .15*meter,'name', 'I','dir','y');

% Gravity
gravity on
gravity y

% Create model
model = MICPModel(G, rock, fluid);

% Boundary condition
f = boundaryFaces(G);
f = f(abs(G.faces.normals(f,1))>eps & (G.faces.centroids(f,1)<-L+2 | ...
                                             G.faces.centroids(f,1)>L-2 ));
fp = G.faces.centroids(f,2)*fluid.rhoWS*norm(gravity);
bc = addBC([], f, 'pressure', fp, 'sat', [0 0]);
bc.o = zeros(size(bc.sat,1), 1);
bc.u = zeros(size(bc.sat,1), 1);
bc.m = zeros(size(bc.sat,1), 1);
bc.b = zeros(size(bc.sat,1), 1);
bc.c = zeros(size(bc.sat,1), 1);

% Well different rates and times
N = 17; % Number of injection changes
M = zeros(N,5); % Matrix where entries per row are: time, rate, o, u, m.
M(1,1) = 15*hour/dt; 
M(1,2) = Q1;
M(2,1) = 26*hour/dt; 
M(2,2) = Q2; 
M(3,1) = 100*hour/dt; 
M(3,2) = Q1;
M(3,3) = fluid.omax;
M(4,1) = 130*hour/dt;
M(4,2) = Q1;
M(5,1) = 135*hour/dt; 
M(5,2) = Q2; 
M(6,1) = 160*hour/dt; 
M(6,2) = Q1;
M(6,4) = fluid.umax;
M(7,1) = 200*hour/dt; 
M(7,2) = Q1;
M(8,1) = 210*hour/dt; 
M(8,2) = Q2;
M(9,1) = 600*hour/dt;
M(9,2) = Q1;
M(9,3) = fluid.omax;
M(10,1) = 630*hour/dt; 
M(10,2) = Q1;
M(11,1) = 650*hour/dt; 
M(11,2) = Q2; 
M(12,1) = 670*hour/dt; 
M(12,2) = Q1;
M(12,4) = fluid.umax;
M(13,1) = 690*hour/dt;
M(13,2) = Q1;
M(14,1) = 710*hour/dt; 
M(14,2) = Q2; 
M(15,1) = 800*hour/dt; 
M(15,2) = Q1;
M(15,4) = fluid.umax;
M(16,1) = 820*hour/dt; 
M(16,2) = Q1;
M(17,1) = 840*hour/dt; 
M(17,2) = Q2; 

% Make schedule
schedule = simpleSchedule(timesteps,'W',W,'bc',bc);
for i=1:N
    schedule.control(i+1)=schedule.control(i);
    schedule.control(i+1).W(1).val=M(i,2)/10;
    schedule.control(i+1).W(2).val=9*M(i,2)/10;
    schedule.control(i+1).W(1).sign=sign(M(i,2));
    schedule.control(i+1).W(2).sign=sign(M(i,2));
    schedule.control(i+1).W(1).o=M(i,3);
    schedule.control(i+1).W(1).u=M(i,4);
    schedule.control(i+1).W(1).m=M(i,5);
    schedule.step.control(M(i,1):end)=i+1;
end    

% Initial condition
state0   = initState(G, W, G.cells.centroids(:,2)*fluid.rhoWS* ...
                                                    norm(gravity), [1, 0]);
state0.o = zeros(G.cells.num,1);
state0.u = zeros(G.cells.num,1);
state0.m = zeros(G.cells.num,1);
state0.b = zeros(G.cells.num,1);
state0.c = zeros(G.cells.num,1);

% Simulate
[~, states] = simulateScheduleADMICP(state0, model, schedule);

% CO2 assesment
statesb=states{600};
statesc=states{800};
statesd=states{nt};
[statese,lr0,~,~]=co2_2Dfls(state0,G,fluid);
[statesf,lr1,~,~]=co2_2Dfls(statesb,G,fluid);
[statesg,lr2,~,~]=co2_2Dfls(statesc,G,fluid);
[statesh,lr3,tt,QCO2]=co2_2Dfls(statesd,G,fluid);

% Write the results to be read in ParaView (GNU Octave)
if exist('OCTAVE_VERSION', 'builtin') ~= 0
    mkdir vtk_micp_2Dfls;
    cd vtk_micp_2Dfls;
    mrsttovtk(G,states,'states','%f');
    %mrsttovtk(G,statesh,'statesh','%f');
    return
end

% Figure 10 paper (MATLAB)
cc(:,1)=[.75:.01:1]';
cc(:,2)=[.75:.01:1]';
cc(:,3)=[.75:-.03:0]';
porosityf=.15-statesb.c-statesb.b;
porosityg=.15-statesc.c-statesc.b;
porosityh=.15-statesd.c-statesd.b;
c=flipud(jet);
c=c(70:1:100,:);
ccc=flipud(jet);
ccc=ccc(70:1:end,:);
set(gcf,'PaperUnits','inches','PaperSize',[9.11 1.85],'PaperPosition', ...
                                                          [0 0 9.11 4.83]);
set(gca,'FontName','Arial');
hold on
n1=subplot(2,4,1);
colormap (n1,cc);
caxis([2e-14 1e-12]);
cb = colorbar; 
title(cb, '$m^2$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.265 .67 .005 .15],'Ticks',[2e-14 1e-12],'FontSize',8);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(a)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
plotCellData(G,K0);
title('Initial permeability','FontSize',8,'FontName','Arial', ...
                                                    'Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax1=axes('position',[.195 .71 .045 .09],'YAxisLocation','right');
box on 
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,K0);
s.EdgeColor = 'none';
colormap (ax1,cc);
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
view(0,270)
n2=subplot(2,4,2);
view(0, 0);
colormap (n2,ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '\%','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.475 .67 .005 .15],'YTick',[0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(b)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,100*(K0-fluid.K(.15-statesb.c-statesb.b))./K0);
s.EdgeColor = 'none';
title('Permeability (phase I MICP)','FontSize',8,'FontName','Arial', ...
                                                    'Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax2=axes('position',[.4 .71 .045 .09],'YAxisLocation','right');
box on 
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,100*(K0-fluid.K(.15-statesb.c-statesb.b))./K0);
s.EdgeColor = 'none';
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
colormap (ax2,ccc);
caxis([0 100]);
view(0,270)
n3=subplot(2,4,3);
view(0, 0);
colormap (n3,ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '\%','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.68 .67 .005 .15],'YTick',[0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(c)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,100*(K0-fluid.K(.15-statesc.c-statesc.b))./K0);
s.EdgeColor = 'none';
title('Permeability (phase II MICP)','FontSize',8,'FontName','Arial', ...
                                                    'Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax3=axes('position',[.605 .71 .045 .09],'YAxisLocation','right');
box on 
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,100*(K0-fluid.K(.15-statesc.c-statesc.b))./K0);
s.EdgeColor = 'none';
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
colormap (ax3,ccc);
caxis([0 100]);
view(0,270)
n4=subplot(2,4,4);
view(0, 0);
colormap (n4,ccc);
caxis([0 100]);
cb = colorbar; 
title(cb, '\%','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.89 .67 .005 .15],'YTick',[0 100]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(d)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,100*(K0-fluid.K(.15-statesd.c-statesd.b))./K0);
s.EdgeColor = 'none';
title('Permeability (phase III MICP)','FontSize',8,'FontName','Arial', ...
                                                    'Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax4=axes('position',[.815 .71 .045 .09],'YAxisLocation','right');
box on 
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,100*(K0-fluid.K(.15-statesd.c-statesd.b))./K0);
s.EdgeColor = 'none';
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
colormap (ax4,ccc);
caxis([0 100]);
view(0,270)
n5=subplot(2,4,5);
view(0, 0);
colormap (n5,c);
caxis([0 75]);
cb = colorbar; 
title(cb,'kg/m$^3$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.265 .2 .005 .15],'YTick',[0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(e)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,fluid.rhoOS*.15*statese.s(:,2));
s.EdgeColor = 'none';
title('CO$_2$ (100 days)','FontSize',8,'FontName','Arial', ...
                                                    'Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax5=axes('position',[.195 .71/3 .045 .09],'YAxisLocation','right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,fluid.rhoOS*.15*statese.s(:,2));
s.EdgeColor = 'none';
colormap (ax5,c);
caxis([0 75]);
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
view(0,270)
n6=subplot(2,4,6);
view(0, 0);
colormap (n6,c);
caxis([0 75]);
cb = colorbar; 
title(cb, 'kg/m$^3$','FontSize',8,'Interpreter','latex','FontName', ...
                                                                  'Arial');
set(cb,'position',[.475 .2 .005 .15],'YTick',[0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(f)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,fluid.rhoOS*porosityf.*statesf.s(:,2));
s.EdgeColor = 'none';
title('CO$_2$ (phase I MICP)','FontSize',8,'FontName', ...
                                            'Arial','Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax6=axes('position',[.4 .71/3 .045 .09],'YAxisLocation','right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,fluid.rhoOS*porosityf.*statesf.s(:,2));
s.EdgeColor = 'none';
colormap (ax6,c);
caxis([0 75]);
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
view(0,270)
n7=subplot(2,4,7);
view(0, 0);
colormap (n7,c);
caxis([0 75]);
cb = colorbar; 
title(cb,'kg/m$^3$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.68 .2 .005 .15],'YTick',[0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(g)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,fluid.rhoOS*porosityg.*statesg.s(:,2));
s.EdgeColor = 'none';
title('CO$_2$ (phase II MICP)','FontSize',8,'FontName', ...
                                            'Arial','Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax7=axes('position',[.605 .71/3 .045 .09],'YAxisLocation','right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,fluid.rhoOS*porosityg.*statesg.s(:,2));
s.EdgeColor = 'none';
colormap (ax7,c);
caxis([0 75]);
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
view(0,270)
n8=subplot(2,4,8);
view(0, 0);
colormap (n8,c);
caxis([0 75]);
cb = colorbar; 
title(cb,'kg/m$^3$','FontSize',8,'Interpreter','latex','FontName','Arial');
set(cb,'position',[.89 .2 .005 .15],'YTick',[0 25 50 75]);
axis([0 L 0 H]);
xlim([0 L])
ylim([0 H])
xlabel({'x [m]'; '(h)'},'FontSize',8,'FontName','Arial');
ylabel('z [m]','FontSize',8,'FontName','Arial');
s=plotCellData(G,fluid.rhoOS*porosityh.*statesh.s(:,2));
s.EdgeColor = 'none';
title('CO$_2$ (phase III MICP)','FontSize',8,'FontName', ...
                                            'Arial','Interpreter','latex');
set(gca,'FontSize',8,'XTick',(0:100:L),'YTick',(0:20:H),'color', ...
                                                'none','FontName','Arial');
view(0, 270);
set(gca,'FontName','Arial');
line([200 100], [60 30], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
line([200 100], [105 130], [0 0],'Color','black','LineStyle','--', ...
                                                            'LineWidth',1);
ax8=axes('position',[.815 .71/3 .045 .09],'YAxisLocation','right');
box on
axis([99.85 100.15 0 H]);
xlim([99.75 100.25])
ylim([30 130])
s=plotCellData(G,fluid.rhoOS*porosityh.*statesh.s(:,2));
s.EdgeColor = 'none';
colormap (ax8,c);
caxis([0 75]);
set(gca,'FontSize',6,'XTick',[99.85 100.15],'YTick',[30 80 130], ...
                                        'color','none','FontName','Arial');
view(0,270)
%print -depsc2 Fig10.eps

% Figures 11a and 11b paper
clear c
clear b
clear u
clear m
clear o
clear K
clear vc
clear v
clear cell_leak
cells =  1:1:G.cells.num;
cell_leak = cells(G.cells.centroids(:,2)<130 & ...
                                    G.cells.centroids(:,2)>30);
for i=1:1:nt
  c(i)=mean(states{i}.c(cell_leak));
  b(i)=mean(states{i}.b(cell_leak));
  m(i)=mean(states{i}.m(cell_leak));
  u(i)=mean(states{i}.u(cell_leak));
  o(i)=mean(states{i}.o(cell_leak));
  Ki=fluid.K(.15-states{i}.c-states{i}.b);
  K(i)=mean(Ki(cell_leak)./K0(cell_leak));
  vc=faceFlux2cellVelocity(G,states{i}.flux(:));
  v(i) = mean(sqrt(sum(vc(cell_leak,:).^ 2, 2)));
end
                              
figure('Units','inches','Position',[0 0 6.83 6.83],...
                                               'PaperPositionMode','auto');
set(gca,'FontName','Arial');
hold on
plot(1:nt,v/max(v),'color',[0 .74 1],'LineWidth',3,'LineStyle','-');
plot(1:nt,m/max(m),'color',[0 .8 0],'LineWidth',3,'LineStyle','-');
plot(1:nt,o/max(o),'color',[1 .5 .9],'LineWidth',3,'LineStyle','-');
plot(1:nt,u/max(u),'color',[1 .9 0],'LineWidth',3,'LineStyle','-');
plot(1:nt,b/max(b),'color',[0 .4 0],'LineWidth',3,'LineStyle',':');
plot(1:nt,c/max(c),'color',[1 .2 .2],'LineWidth',3,'LineStyle',':');
plot(1:nt,K,'color',[0 0 0],'LineWidth',3,'LineStyle',':');
line([600 600], [0 1], [0 0],'Color',[0 0 0],'LineStyle','--', ...
                                                           'LineWidth',1);
line([800 800], [0 1], [0 0],'Color',[0 0 0],'LineStyle','--', ...
                                                           'LineWidth',1);                                                       
hold off

text(250,1.02,'Phase I','FontSize',11,'Interpreter','latex', ...
                                                       'FontName','Arial');
text(650,1.02,'Phase II','FontSize',11,'Interpreter','latex', ...
                                                       'FontName','Arial');
text(825,1.02,'Phase III','FontSize',11,'Interpreter','latex', ...
                                                       'FontName','Arial');                                                   
xlim([0 nt]);
xlabel({'Time [h]';'(a)'},'FontSize',11,'Interpreter','latex');        
ylabel('[$-$]','FontSize',11,'Interpreter','latex');
h=legend('$v_w/0.0070\textrm{ m/s}$','$c_m/0.0019\textrm{ kg/m}^3$',...
'$c_o/0.0085\textrm{ kg/m}^3$','$c_u/141 \textrm{ kg/m}^3$',...
'$\phi_b/0.0004$','$\phi_c/0.0363$','$K/10^{-12}\textrm{ m}^2$',...
                    'Interpreter','latex','FontSize',11);
rect = [0.37, 0.35, .2, .25];
set(h, 'Position', rect);               
set(gca,'FontSize',11,'FontName','Arial','XTick',[0:100:1000],'YGrid', ...
                                                      'on', 'XGrid', 'on');
%print -depsc2 Fig11a.eps

figure('Units','inches','Position',[0 0 6.83 6.83], ...
                                               'PaperPositionMode','auto');
set(gca,'FontName','Arial');
hold on
plot(tt(1,:),100*lr0/QCO2,'color',[1 .2 .2],'LineWidth',9,'LineStyle','-');
plot(tt(1,:),100*lr1/QCO2,'color',[1 .5 0],'LineWidth',9,'LineStyle','-');
plot(tt(1,:),100*lr2/QCO2,'color',[0.61 0.61 0.61],'LineWidth',9, ...
                                                          'LineStyle','-');
plot(tt(1,:),100*lr3/QCO2,'color',[0 0 0],'LineWidth',9,'LineStyle','-');
hold off
xlim([0 100]);
ylim([0 60]);
xlabel({'Time [d]';'(b)'},'FontSize',11,'Interpreter','latex');        
ylabel('CO$_2$ leakage rate/injection rate [\%]','FontSize',11, ...
                                                    'Interpreter','latex');
grid on
legend('Without MICP','Phase I MICP','Phase II MICP','Phase III MICP', ...
                                                        'Location','best');
set(gca,'FontSize',11,'FontName','Arial','XTick',(0:20:100), ...
                                                        'YTick',(0:10:60));
%print -depsc2 Fig11b.eps