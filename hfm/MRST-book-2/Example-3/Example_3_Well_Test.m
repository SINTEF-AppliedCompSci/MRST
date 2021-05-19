%% LOAD MRST MODULES
mrstModule add hfm ad-blackoil ad-core ad-props mrst-gui

%% LOAD APODI 2 DATA
load Apodi2.mat

%% SETUP GRID
celldim=[220 220 1];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);
km=1*milli*darcy;
G.rock=makeRock(G,km,0.25);

%% EDFM PRE-PROCESSING
tol=1e-6;
[G,fracplanes]=EDFMgrid(G,fracplanes,'Tolerance',tol);
G=fracturematrixNNC3D(G,tol);
[G,fracplanes]=fracturefractureNNCs3D(G,fracplanes,tol);

%% SETUP FLUID MODEL WITH OIL PROPERTIES
pRef = 100*barsa;
fluid = initSimpleADIFluid('phases','W',...
'mu' , 8* centi*poise     , ...
'rho', 700* kilogram/meter^3, ...
'pRef', pRef, ...
'c', 1e-5/barsa);

%% SETUP WATER MODEL
gravity reset off
model = WaterModel(G, G.rock, fluid);
model.operators = setupEDFMOperatorsTPFA(G, G.rock, tol);

%% PLOT GRID TO HELP DETERMINE APPROXIMATE WELL LOCATION
[hm,hf]=plotEDFMgrid(G);
hm.EdgeAlpha=0;
hold on
scatter(27.16, 121, 10, 'r', 'filled');
view(0,90);
box on
xlabel('x [m]')
ylabel('y [m]')

%% LOCATE NEAREST FRACTURE CELL
fraccellcent = [27.16,121,0.5]; % Approximate location of well
fraccell=find(all(abs(G.cells.centroids-fraccellcent)<0.5,2));
fraccell=fraccell(fraccell>G.Matrix.cells.num);
assert(length(fraccell)==1,'more/less than one fraccell detected');

%% SET UP SINK
src=[];
src=addSource(src,fraccell,-1*(meter^3)/day);
src.sat=ones(size(src.cell));

%% RUN SIMULATION
state  = initResSol(G, pRef);
dt = diff(10.^(-8.5+(1:60)*0.25))'; % log spaced timesteps
schedule = simpleSchedule(dt,'src',src);
[~, states] = simulateScheduleAD(state, model, schedule);

%% WELL TEST DIAGNOSTIC PLOT
cellnum=src.cell;
pvalues=zeros(length(states),1);

for i=1:length(states)
	pvalues(i)=states{i}.pressure(cellnum);   
end

tottime=cumsum(dt);

% calc derivative
dp=-diff(pvalues);
dlnt=diff(log(tottime));
dp_dlnt=dp./dlnt;

% get times
midtime=0.5*(tottime(1:end-1)+tottime(2:end));

loglog(midtime,dp_dlnt);
grid on
xlabel('Time (s)')
ylabel('Pressure derivative (Pa)')

%% PLOT MATRIX RESULTS
figure;
plotToolbar(G,states,'lockCaxis',true, 'field', 'pressure');
caxis([50 100]*barsa);
colorbar('EastOutside');
[hm,hf] = plotEDFMgrid(G);
hm.EdgeAlpha=0;
hf.EdgeAlpha=0.25;
xlabel('x [m]')
ylabel('y [m]')
view(0, 90);
axis equal tight
box on

%% PLOT FRACTURE RESULTS
figure;
Gplot = createGplot(G,1,fracplanes);
plotToolbar(Gplot,states,1+G.Matrix.cells.num:G.cells.num,...
    'lockCaxis',true, 'field', 'pressure');
caxis([50 100]*barsa);
colorbar('EastOutside');
xlabel('x [m]')
ylabel('y [m]')
view(0, 90);
axis equal tight
box on
