%% LOAD MRST MODULES
mrstModule add hfm ad-blackoil ad-core ad-props

%% FRACTURE NETWORK STATISTICAL PARAMETERS
physdim=[100 100 100]; % domain size

% Fracture sets 1 (red), 2(green), 3(yellow)
fracinput1=struct('vertices',10,...
'P32',0.025*(1/meter),...
'normal',struct('direction',[1 0 0],'K',10),...
'size',struct('minsize',5*meter,'maxsize',20*meter,'exponent',1.5),...
'perm',1e7*milli*darcy,...
'poro',1,...
'aperture',0.001*meter,...
'circle',true);
fracinput2=fracinput1; fracinput2.normal.direction=[0 0 1];
fracinput3=fracinput1; fracinput3.normal.direction=[0 1 0];

% Exclusion zone is cylindrical with radius being (1+exclzonemult) times
% the fracture radius; height is exclzonemult times the fracture radius.
exclzonemult=0.01;

%% DFN GENERATION
tol=10^-5; % tolerance for comparison of doubles
fracplanes=[]; % empty list of fracture objects

% set 1
[fracplanes,~]=DFNgenerator(fracplanes,fracinput1,physdim,...
exclzonemult,tol);
% set 2
[fracplanes,~]=DFNgenerator(fracplanes,fracinput2,physdim,...
exclzonemult,tol);
% set 3
[fracplanes,~]=DFNgenerator(fracplanes,fracinput3,physdim,...
exclzonemult,tol);

%% PLOT DFN
colourchoice=['r','g','y'];
figure;
hold on;
plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
axis equal tight
view(45,30);
for i=1:3
    index = find(vertcat(fracplanes.SetID)==i);
    C=colourchoice(i);
    for j=index'
        X=fracplanes(j).points(:,1);
        Y=fracplanes(j).points(:,2);
        Z=fracplanes(j).points(:,3);
        fill3(X,Y,Z,C);
    end
    xlabel('x','Interpreter','latex')
    ylabel('y','Interpreter','latex')
    zlabel('z','Interpreter','latex')
end

%% SETUP GRID
celldim=[25 25 25];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
km=10*milli*darcy;
G.rock=makeRock(G,km,0.25);

%% EDFM PRE-PROCESSING
tol=1e-6;
[G,fracplanes]=EDFMgrid(G,fracplanes,'Tolerance',tol);
G=fracturematrixNNC3D(G,tol);
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

%% SETUP FLUID MODEL WITH WATER PROPERTIES
pRef = 100*barsa;
fluid = initSimpleADIFluid('phases','W',...
'mu' , 1* centi*poise     , ...
'rho', 1000* kilogram/meter^3, ...
'pRef',    pRef, ...
'c',       0/barsa);

%% SETUP WATER MODEL
gravity reset off
model = WaterModel(G, G.rock, fluid);
model.operators = setupEDFMOperatorsTPFA(G, G.rock, tol);
model.stepFunctionIsLinear=true;

%% INITIAL AND BOUNDARY CONDITIONS
% Initial condition
state  = initResSol(G, pRef); 

% Find fracture cell faces at domain boundary
boundfaces=findfracboundaryfaces(G,tol);

% Set pressure differential on opposing boundaries in the x-direction
deltaP = 50*barsa;
bc=[];
bc=pside(bc,G.Matrix,'East',pRef);
matwestfaces=bc.face;
bc=pside(bc,G.Matrix,'West',pRef + deltaP);
bc=addBC(bc,boundfaces.East,'pressure',100*barsa);
bc=addBC(bc,boundfaces.West,'pressure',150*barsa);
bc.sat=ones(size(bc.face));

%% SETUP SCHEDULE
schedule = simpleSchedule(1,'bc',bc);

%% SIMULATE
[~, states,~] = simulateScheduleAD(state, model, schedule);

%% PLOT RESULTS
figure;
plotCellData(G,states{1}.pressure,1:G.Matrix.cells.num,...
    'FaceAlpha',0.5,'EdgeAlpha',0.1);
plotCellData(G,states{1}.pressure,G.Matrix.cells.num+1:G.cells.num);
view(30,45);
caxis([100 150]*barsa);
colorbar('EastOutside');
axis equal tight
box on

%% CALCULATE EQUIVALENT PERMEABILITY
% Determine flux through western boundary
westfaces=[matwestfaces;boundfaces.West'];
waterflux=sum(abs(states{end}.flux(westfaces,1)));

% Darcy's law
k_eq=waterflux*fluid.muW(1)*physdim(1)/(physdim(1)*physdim(2)*deltaP);
disp(convertTo(k_eq,milli*darcy));