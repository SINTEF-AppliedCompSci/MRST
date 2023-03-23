%% This will generate Figures 10 to 13 in SPE-201243-PA when the fractures 
%  are highly conductive
%% PEDFM Validation case of + fracture of single phase flow
% This code is used to generate the plots in Figures 10 to 13 of 
% Olorode et al., 2020, which are based on Tene et al., 2017 Figures 3 & 4
clear 
clc
close all
opt = struct('nkr',       2, ...
    'shouldPlot', 1 );
%     opt = merge_options(opt, varargin{:});
%% Load necessary modules
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add ad-blackoil;     % AD blackoil solver
mrstModule add ad-core;         % AD core module
mrstModule add ad-props;        % AD properties

%% Basic Parameters
tol=1e-5;
physdim = [9 9 3]*meter;
griddim = [27 27 3];
cellsize = physdim./griddim;

%% Create Box Reservoir Grid
G_matrix = tensorGrid(0:cellsize(1):physdim(1),...
                      0:cellsize(2):physdim(2),...
                      0:cellsize(3):physdim(3));
G_matrix = computeGeometry(G_matrix);
G_matrix.rock=makeShaleRock(G_matrix,1*darcy,0.3);

%% Create Fracture System
fracplanes = struct;
fracplanes(1).points=...
    [2,4.5,0;
     7,4.5,0;
     7,4.5,3;
     2,4.5,3];

fracplanes(2).points=...
   [4.5,2,0;
    4.5,7,0;
    4.5,7,3;
    4.5,2,3];
  
for i=1:numel(fracplanes)
     fracplanes(i).aperture = 1/25; %1*micro*meter;
     fracplanes(i).poro = 0.3;
     fracplanes(i).perm = 1e3*darcy;
%      fracplanes(i).perm = 1e-8*darcy;
end


 %% Create Wells
wells = struct;

%% visualize to check before pre-process
checkIfCoplanar(fracplanes)

G=G_matrix;

%% EDFM PreProcessing
[G,fracplanes]=EDFMshalegrid(G,fracplanes,...
        'Tolerance',tol,'plotgrid',false,'fracturelist',1:numel(fracplanes));
  %-Frac-Matrix NNCs
G = fracturematrixShaleNNC3D(G,tol);
  %-Frac-Frac NNCs
[G,fracplanes]=fracturefractureShaleNNCs3D(G,fracplanes,tol);
  %-Well-Fracs NNCs
[G,wells] = wellfractureNNCs3D(G,fracplanes,wells,tol);

  % OMO: Projection-based NNCs % Switch from EDFM to pEDFM
Ge = G; %EDFM
G = pMatFracNNCs3D(G,tol); %pEDFM

figure,
plotFracSystem(G_matrix,fracplanes,wells);

%% Define resolved explict grid
griddim_r = [225 225 2];
cellsize_r = physdim./griddim_r;

Gr = tensorGrid(0:cellsize_r(1):physdim(1),...
                      0:cellsize_r(2):physdim(2),...
                      0:cellsize_r(3):physdim(3));
Gr = computeGeometry(Gr);

Gr.rock=makeShaleRock(Gr,1*darcy,0.3);


%Set fracture perm
IJ_hFrac=ones(175-51+1,2)*113; IJ_hFrac(:,1)=51:175;
IJ_vFrac=ones(175-51+1,2)*113; IJ_vFrac(:,2)=51:175;
IJ_Fracs=[IJ_hFrac;IJ_vFrac];
FracCellIDs_2D = sub2ind(Gr.cartDims(1:2), IJ_Fracs(:,1), IJ_Fracs(:,2));

fractureCells=[];
NumGridLayer=griddim_r(1)*griddim_r(2);
for layerID = 1:griddim_r(3)
    fractureCells=[fractureCells; FracCellIDs_2D+(layerID-1)*NumGridLayer];
end
Gr.rock.perm(fractureCells) = fracplanes(1).perm; %Frac perm

figure,
show = false([Gr.cells.num, 1]);
show(fractureCells) = true;% Hide well cell
clf; 
plotCellData(Gr,convertTo(Gr.rock.perm,darcy),show,'EdgeColor', 'None','facealpha',0.5);
plotCellData(Gr,convertTo(Gr.rock.perm,darcy),~show,'EdgeColor', 'None','facealpha',0.0);
view(3); colorbar; axis equal tight

    
%% Define fluid properties
% Define a single-phase fluid model without capillarity.

fluid = initSimpleADIFluid('phases','W',       ... % Fluid phase: water
                           'mu',  1*centi*poise, ... % Viscosity
                           'rho', 1000,          ... % Surface density [kg/m^3]
                           'c',   1e-4/barsa,    ... % Fluid compressibility
                           'cR',  1e-5/barsa     ... % Rock compressibility
                           );

%% Define single-phase compressible flow model with pEDFM operators
% We define a three-phase black-oil model without dissolved gas or vaporized
% oil. This is done by first instantiating the blackoil model, and then
% manually passing in the internal transmissibilities and the topological
% neighborship from the embedded fracture grid.

gravity reset off
model = WaterModel(G,[], fluid);
model.operators = setupPEDFMOpsTPFA(G, G.rock, tol);

model_e = WaterModel(Ge,[], fluid);
model_e.operators = setupShaleEDFMOpsTPFA(Ge, Ge.rock, tol);

model_r = WaterModel(Gr,[], fluid);
model_r.operators = setupOperatorsTPFA(Gr, Gr.rock);

s0 = (1);
bc = [];
bc  = pside(bc, G, 'LEFT', 10*barsa,'sat', s0);
bc  = pside(bc, G, 'RIGHT', 0*barsa,'sat', s0);

bce = [];
bce  = pside(bce, Ge, 'LEFT', 10*barsa,'sat', s0);
bce  = pside(bce, Ge, 'RIGHT', 0*barsa,'sat', s0);

bcr = [];
bcr  = pside(bcr, Gr, 'LEFT', 10*barsa,'sat', s0);
bcr  = pside(bcr, Gr, 'RIGHT', 0*barsa,'sat', s0);

clf, figure
plotGrid(G_matrix,'FaceColor', 'none'); view(3);
plotFaces(G_matrix, bc.face(strcmp(bc.type,'pressure')), 'r','facealpha',0.7);
plotFaces(G_matrix, bc.face(strcmp(bc.type,'pressure')), 'r','facealpha',0.7);
hold on;

hf={};
for i = 1:length(fracplanes)
    X=fracplanes(i).points(:,1);
    Y=fracplanes(i).points(:,2);
    Z=fracplanes(i).points(:,3);
    hf{end+1}=fill3(X,Y,Z,'b');
    set(hf{end},'facealpha',0.7);
end
set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Set up initial state and schedule
% We set up a initial state with the reference pressure and a mixture of
% water and oil initially. We also set up a simple-time step strategy that
% ramps up gradually towards 30 day time-steps.

s0 = (1);
pRef=5*barsa;
state  = initResSol(G, pRef, s0);
state_e = initResSol(Ge, pRef, s0);
state_r = initResSol(Gr, pRef, s0);

totTime = 60*day;
nSteps = 10;
dt = rampupTimesteps(totTime, 30*day, nSteps);

schedule = simpleSchedule(dt, 'bc', bc);
schedule_e = simpleSchedule(dt, 'bc', bce);
schedule_r = simpleSchedule(dt, 'bc', bcr);

%% Simulate problem
[ws, states, reports] = simulateScheduleAD(state, model, schedule, ...
   'afterStepFn', getPlotAfterStep(state, model, schedule));

[ws_e, states_e, reports_e] = simulateScheduleAD(state_e, model_e, schedule_e, ...
   'afterStepFn', getPlotAfterStep(state_e, model_e, schedule_e));

[ws_r, states_r, reports_r] = simulateScheduleAD(state_r, model_r, schedule_r, ...
   'afterStepFn', getPlotAfterStep(state_r, model_r, schedule_r));

%% plotting
figure, 
plotToolbar(G, states)
view(40,30);
axis tight equal;

%% Pressrue plot
figure,
plotCellData(G, convertTo(states{end,1}.pressure, barsa()), 'EdgeColor', 'k');
xlabel('x'), ylabel('y'), zlabel('Depth');
view(3); axis tight;colormap('jet');
h=colorbar; caxis([0 10]); 
h.Label.String = 'Pressure (bar)';
set(gca, 'CameraPosition', [-24.1107  -33.7559  -19.0289]);


%% Pressure contour comparsion

% Figure and figure settings
figure('Position',[200 460 1200 400]);

[X,Y]=meshgrid(cellsize_r(1)/2:cellsize_r(1):physdim(1),cellsize_r(2)/2:cellsize_r(2):physdim(2));
nLayer=1;
NumGridLayer=griddim_r(1)*griddim_r(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;

pres_r=reshape(states_r{end,1}.pressure(startId:endId),[griddim_r(1),griddim_r(2)]);

subplot(1,3,1);
surf(X,Y,convertTo(pres_r, barsa()),'EdgeAlpha',0.1), view(130,20)
colormap('jet');axis tight;
title('Explicit solution pressure')

[X,Y]=meshgrid(cellsize(1)/2:cellsize(1):physdim(1),cellsize(2)/2:cellsize(2):physdim(2));
nLayer=1;
NumGridLayer=griddim(1)*griddim(2);
startId=(nLayer-1)*NumGridLayer+1;
endId=(nLayer)*NumGridLayer;

pres=reshape(states{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,3,2);
surf(X,Y,convertTo(pres, barsa()),'EdgeAlpha',0.3), view(130,20)
colormap('jet');axis tight;
title('pEDFM pressure')

pres_e=reshape(states_e{end,1}.pressure(startId:endId),[griddim(1),griddim(2)]);

subplot(1,3,3);
surf(X,Y,convertTo(pres_e, barsa()),'EdgeAlpha',0.3), view(130,20)
colormap('jet');axis tight;
title('EDFM pressure')




%% Pressure overline comparsion
IJ_hFrac=ones(griddim(1),2)*((griddim(1)-1)/2+1); IJ_hFrac(:,1)=1:griddim(1);
CenterCellIDs = sub2ind(G.cartDims(1:2), IJ_hFrac(:,1), IJ_hFrac(:,2));

IJ_hFrac_r=ones(225,2)*113; IJ_hFrac_r(:,1)=1:225;
CenterCellIDs_r = sub2ind(Gr.cartDims(1:2), IJ_hFrac_r(:,1), IJ_hFrac_r(:,2));

x_line=IJ_hFrac(:,1)*cellsize(1)-cellsize(1)/2;

pres_line=states{end,1}.pressure(CenterCellIDs);
pres_line_e=states_e{end,1}.pressure(CenterCellIDs);


x_line_r=IJ_hFrac_r(:,1)*cellsize_r(1)-cellsize_r(1)/2;
pres_line_r=states_r{end,1}.pressure(CenterCellIDs_r);


figure('rend','painters','pos',[10 10 640 500]);
title('Pressure solution over the horizontal center line');
plot(x_line_r, convertTo(pres_line_r, barsa()),'k-', 'LineWidth', 2,'DisplayName','Explicit');
hold on;
plot(x_line, convertTo(pres_line, barsa()),'o','Color','r', ...
    'LineWidth', 2,'DisplayName','pEDFM');
plot(x_line, convertTo(pres_line_e, barsa()),'^','Color','b', ...
    'LineWidth', 2,'DisplayName','EDFM');
hold off;

xlim([-1e-3 9]);
ylim([-1e-3 10]);

set(gca,'FontSize',15);
xlabel('X, m');
ylabel('Pressure, Bar');
legend; axis tight;