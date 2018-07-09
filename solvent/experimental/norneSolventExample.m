mrstModule add ad-core ad-eor ad-blackoil ad-props blackoil-sequential matlab_bgl

gravity reset on

%%

mrstModule add deckformat
grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);
tmp = load(fullfile(mrstDataDirectory(), 'norne.mat'));
G = tmp.G;
start = 4;
stop = 22;
[ii, jj, kk] = gridLogicalIndices(G);
G = extractSubgrid(G, find(kk >= start & kk <= stop));

G = computeGeometry(G);
rock = grdecl2Rock(grdecl, G.cells.indexMap);

% W = tmp.W;
% for i = 1:numel(W)
%     rem = W(i).cells > G.cells.num;
%     W(i).cells(rem) = [];
%     W(i).r(rem) = [];
%     W(i).dir(rem) = [];
%     W(i).WI(rem) = [];
%     W(i).dZ(rem) = [];
%     W(i).status = true;
%     W(i).cstatus = true;
%     W(i).lims = [];
%     if W(i).sign < 0
%         W(i).type = 'bhp';
%         W(i).val = 0;
%     end
% end

%%


inj = [9, 15; ...
        26, 15; ...
       36, 80; ...
       10, 85; ...
       24, 30; ...
       14, 52; ...
       18, 80;
       23, 66];

nInj = size(inj,1);
W = [];
T = 5*year;
pv = 0.5*sum(poreVolume(G, rock));
rate = (pv/T)/nInj;
for i = 1:nInj
    W = verticalWell(W, G, rock, inj(i, 1), inj(i, 2), [], ...
        'comp_i', [0,0,0,1],...
        'type', 'rate', 'val', rate);
end
   
prod = [10, 66; ...
        12, 32; ...
        22, 49; ...
        13, 91; ...
        37, 95;
        35, 64];
nProd = size(prod,1);
for i = 1:nProd
    W = verticalWell(W, G, rock, prod(i, 1), prod(i, 2), [], ...
        'comp_i', [1,0,0,0], ...
        'type', 'bhp', 'val', 80*barsa);
end
    
%%

gravity reset on

fluid = initSimpleADIFluid('n'     , [2, 2, 2], ...
                           'rho'   , [1000, 800, 100]*kilogram/meter^3, ...
                           'phases', 'WOG', ...
                           'mu'    , [1, 10, 2]*centi*poise);

sOres_i = 0.38;
sOres_m = 0.08;
fluid = addSolventProperties(fluid, 'n', 2, ...
                                    'rho', 100*kilogram/meter^3, ...
                                    'mixPar', 2/3, ...
                                    'mu'    , 1*centi*poise, ...
                                    'sOres_i', sOres_i, ...
                                    'sOres_m', sOres_m);
                                
model = FourPhaseSolventModel(G, rock, fluid);
model.extraStateOutput = true;

nStep = 200;

dT = repmat(T/nStep, nStep,1);

schedule = simpleSchedule(dT, 'W', W);
sg = 0.1;
sg = 0;
state0 = initResSol(G, 100*barsa, [1-sOres_i-sg sOres_i sg 0]);
state0.wellSol = initWellSolAD(W, model, state0);

% nls = NonLinearSolver('useLineSearch', true);
nls = NonLinearSolver('useLineSearch', false);

%%

% getHandler = @(name) ResultHandler('dataFolder', ['norneWAG', ...
%     '_nStep', num2str(nStep), '_nCycles', num2str(nCycles), ...
%     '_sOres_i', num2str(sOres_i), '_sOres_m', num2str(sOres_m)], ...
%     'dataPrefix', [name, '_step'], 'cleardir', false);
% 
% stateHandler = getHandler('state');
% repHandler = getHandler('rep');

%%

[ws, states, reports] = simulateScheduleAD(state0, model, schedule, 'nonlinearSolver', nls);

%%

plotToolbar(G, states);
