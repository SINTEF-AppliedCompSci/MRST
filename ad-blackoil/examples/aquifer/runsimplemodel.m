mrstModule add mrst-gui diagnostics deckformat opm_gridprocessing coarsegrid ...
    ad-props ad-core ad-blackoil 
mrstVerbose true
gravity on

% fn = ...
dirname = mfilename('fullpath');
dirname = fileparts(dirname);

fn = fullfile(dirname, '2D_OILWATER_AQ_ECLIPSE_MOD.DATA');
deck = readEclipseDeck(fn);
deck = convertDeckUnits(deck);

% G = initEclipseGrid(deck, 'SplitDisconnected', false);
G = initEclipseGrid(deck);
G = computeGeometry(G(1));

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
% rock.perm = 1e9*rock.perm;

% end-point scaling, so include grid
fluid = initDeckADIFluid(deck, 'G', G);

output = processAquifer(deck, G);
aquifers     = output.aquifers;
aquind       = output.aquind; 
initval      = output.initval;
aquiferprops = output.aquiferprops;

simplemodel = AquiferBlackOilModel(G, rock, fluid, aquifers, aquind, aquiferprops, ...
                                   'modeltype', 'oilwater');
simplebomodel = TwoPhaseOilWaterModel(G, rock, fluid);
if isfield(deck.RUNSPEC, 'VAPOIL')
    simplemodel.vapoil = deck.RUNSPEC.VAPOIL;
    simplebomodel.vapoil = deck.RUNSPEC.VAPOIL;
else
    simplemodel.vapoil = 0;        
    simplebomodel.vapoil = 0;        
end
if isfield(deck.RUNSPEC, 'DISGAS')
    simplemodel.disgas = deck.RUNSPEC.DISGAS;
    simplebomodel.disgas = deck.RUNSPEC.DISGAS;
else
    simplemodel.disgas = 0;        
    simplebomodel.disgas = 0;        
end
simplemodel.FacilityModel = selectFacilityFromDeck(deck, simplemodel);
simplebomodel.FacilityModel = selectFacilityFromDeck(deck, simplebomodel);

schedule = convertDeckScheduleToMRST(simplemodel, deck);
schedulebo = convertDeckScheduleToMRST(simplebomodel, deck);

initState = initStateDeck(simplemodel, deck);
initState.aquiferpressures = initval.pressures;
initState.aquifervolumes   = initval.volumes;

initStatebo = initStateDeck(simplebomodel, deck);

[wellSols, states, schedulereport] = simulateScheduleAD(initState, simplemodel, ...
                                                  schedule);
