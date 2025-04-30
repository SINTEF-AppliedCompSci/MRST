clear;
close all;
set(0,'defaultfigurecolor','w')
mrstModule add ad-core ad-props compositional ecpa dual-porosity dual-porosity-permeability deckformat mrst-gui
gravity on

%% Choose between the dual porosity or dual porosity-dual permeability models
% model_name = 'DP';
model_name = 'DPDP';

%% Set up grid
celldim = [40, 19, 10];
physdim = [2000 950 100];
nx = 0:50:physdim(1);
ny = 0:50:physdim(2);
nz = 0:10:physdim(3);

G = tensorGrid(nx, ny, nz);
G = computeGeometry(G);

figure;
plotGrid(G,'facealpha',0,'edgealpha',0.5);
view(30,-40);axis tight
xlabel('Length (m)','FontSize',12);
ylabel('Width (m)','FontSize',12);
zlabel('Height (m)','FontSize',12);

%% Set up rock properties
MatrixPerm = 0.1*milli*darcy;
FracturePerm = 5*milli*darcy;

rock_matrix = makeRock(G, [MatrixPerm,MatrixPerm,0.1*MatrixPerm], 0.05);
rock = makeRock(G, FracturePerm, 0.01);

%% Set up fluid
mixture = ECPATableCompositionalMixture({'CarbonDioxide','Methane'});
T = 273.15 + 85.74; p = 5.5E6;  z = [0 1];     % p, T, z
bic = eCPAreadBinaryInteraction(mixture,T(1));
ap = eCPAreadAssociationParameter(mixture, T(1));
mixture = setBinaryInteraction(mixture, bic);
mixture = setAssociationParameter(mixture, ap);

eCPA = ECPAEquationOfStateModel([], mixture, 'eCPA');
[~, ~, ~,~, ~,rhoOSC,rhoGSC] = eCPAstandaloneFlash(barsa, 293.15, [1 0], eCPA);
fluid_matrix = initSimpleADIFluid('phases', 'OG', 'blackoil', false, 'rho', [rhoOSC,rhoGSC],'n', [2, 2]);
fluid_fracture = initSimpleADIFluid('phases', 'OG', 'blackoil', false, 'rho', [rhoOSC,rhoGSC],'n', [2, 2]);
fluid_matrix.relPermScal = 1;
fluid_fracture.relPermScal = 1;
%% Set the DP model. Here, a single-phase model is used. Rock and fluid
% are sent in the constructor as a cell array {fracture,matrix}
diagonal_backend = DiagonalAutoDiffBackend('modifyOperators', false, 'rowMajor', false, 'useMex', false);
if strcmp(model_name, 'DP')
    natural= ECPANaturalVariablesCompositionalModelDP(G, {rock,rock_matrix},...
        {fluid_fracture,fluid_matrix},...
        mixture,...
        'water',false, 'oil',true, 'gas', true,...
        'liquidPhase', 'O', 'vaporPhase', 'G',...
        'AutoDiffBackend', diagonal_backend);
    overall = ECPAOverallCompositionCompositionalModelDP(G, {rock,rock_matrix},...
        {fluid_fracture,fluid_matrix},...
        mixture,...
        'water',false, 'oil',true, 'gas', true,...
        'liquidPhase', 'O', 'vaporPhase', 'G',...
        'AutoDiffBackend', diagonal_backend);
elseif strcmp(model_name, 'DPDP')
    natural= ECPANaturalVariablesCompositionalModelDPDP(G, {rock,rock_matrix},...
        {fluid_fracture,fluid_matrix},...
        mixture,...
        'water',false, 'oil',true, 'gas', true,...
        'liquidPhase', 'O', 'vaporPhase', 'G',...
        'AutoDiffBackend', diagonal_backend);
    overall = ECPAOverallCompositionCompositionalModelDPDP(G, {rock,rock_matrix},...
        {fluid_fracture,fluid_matrix},...
        mixture,...
        'water',false, 'oil',true, 'gas', true,...
        'liquidPhase', 'O', 'vaporPhase', 'G',...
        'AutoDiffBackend', diagonal_backend);
end
natural = imposeRelpermScaling(natural, 'SOGCR', 0);
overall = imposeRelpermScaling(overall, 'SOGCR', 0);
overall = overall.validateModel();
natural = natural.validateModel();

%% Setting transfer function. This step is important to ensure that fluids
% will be transferred from fracture to matrix (and vice-versa). There are 
% several options of shape factor (see folder
% transfer_models/shape_factor_models/) that could be set in this
% constructor. The second argument is the matrix block size. Another
% possible transfer function to be used in this model would be:
%       model.transfer_model_object = EclipseTransferFunction();

overall.transfer_model_object = KazemiTwoPhaseTransferFunction('KazemiShapeFactor',...
                                                                 [10,10,10]);
natural.transfer_model_object = KazemiTwoPhaseTransferFunction('KazemiShapeFactor',...
                                                                 [10,10,10]);
%% Initializing state
[L1, ~, ~,Z_L, Z_V] = eCPAstandaloneFlash(p, T, z, eCPA);
sV = Z_V.*(1-L1)./(Z_V.*(1-L1)+Z_L.*L1);
s0 = [1-sV, sV];

[~, ~, ~,~, ~,rho0] = eCPAstandaloneFlash(p, T, z, eCPA);  
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rho0, [z_0, z_max], p);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
state = eCPAinitCompositionalState(overall, p0, T, s0, z, true);
state = overall.validateState(state);
initGuess = state;
initGuess = overall.computeFlash(initGuess);

% prod well
wellprod = celldim(1)*(celldim(2)+1)/2;
cellprod = (wellprod:celldim(1)*celldim(2):celldim(3)*celldim(1)*celldim(2));
W = addWell([], G, rock, cellprod, ...
    'comp_i', [0, 1],'Name', ['Prod',num2str(1)], 'Val',...
    5.5e6, 'sign', -1, 'Type', 'bhp','Radius', 0.1,'Dir','z');
W(1).components = z;
% inj well
wellinj = celldim(1)*(celldim(2)-1)/2+1;
cellinj = (wellinj:celldim(1)*celldim(2):celldim(3)*celldim(1)*celldim(2));
W = addWell(W, G, rock, cellinj, ...
    'comp_i', [0, 1],'Name', ['Inj',num2str(2)], 'Val',...
    20*1e4/day, 'sign', 1, 'Type', 'rate','Radius', 0.1,'Dir','z');
W(2).components = [1 0];
plotWell(G, W, 'height', 10)

t1 = 30*year();
dt = rampupTimesteps(t1,90*day(), 12);
schedule = simpleSchedule(dt, 'W', W);

% Simulate problem
nls = NonLinearSolver();
[ws, states, reports] = simulateScheduleAD(state, natural, schedule, 'nonlinearsolver', nls);
% [wso, stateo, reporto] = simulateScheduleAD(state, overall, schedule, 'nonlinearsolver', nls);

figure, 
plotToolbar(G, states)
view(0,0);
axis tight;
colorbar


res=[];
t=schedule.step.val;
for i = 1:numel(ws)
    mCO2 = ws{i}(1).CarbonDioxide;
    mCH4 = ws{i}(1).Methane;
    res=[res;mCH4,mCO2];
end
res=abs(res).*24.*3600;

timeSim = cumsum(schedule.step.val)./year();

CMG = [2.73973E-05	182.5377655	0
0.486673803	2964.670654	0
1.027456526	13851.95606	0
2.073533729	33382.42188	4362.305176
2.963633162	23196.27148	62345.71484
3.949934532	8522.442383	145820.0469
4.936235901	3373.783447	201296.25
6.169112614	2262.033447	245287.8906
7.155413984	2534.8396	269873
8.141715353	3054.803955	288530.625
9.128016723	3569.528076	302973
10.11431809	3988.707275	314297.0625
11.10061946	4307.914063	323212.3125
12.08692083	4545.197266	330231.8125
13.0732222	4724.123535	335741.2813
14.05952357	4861.743164	340051.4375
15.04582494	4969.604981	343413.9063
16.03212631	5052.404297	346040.8438
17.01842768	5113.910156	348098.0313
18.00472905	5157.04541	349715.7188
18.99103042	5184.850098	350993.4063
19.97733179	5199.80127	352008.6875
22.19650987	5197.67041	353621.5938
24.16911261	5168.121094	354549.5313
26.14171535	5123.51416	355195.5938
28.11431809	5070.286621	355669.6875
30.02191781	5004.671387	356097.75];

figure
plot(timeSim,res(:,1));
hold on 
plot (CMG(:,1),CMG(:,2),'ro');
legend('MRST-CPA', 'CMG')
figure
plot(timeSim,res(:,2));
hold on 
plot (CMG(:,1),CMG(:,3),'ro');
legend('MRST-CPA', 'CMG')
