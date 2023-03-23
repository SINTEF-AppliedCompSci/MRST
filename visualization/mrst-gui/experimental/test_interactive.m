if 0
    grdecl = fullfile(ROOTDIR, 'examples', 'data', 'SAIGUP', 'SAIGUP.GRDECL');
    grdecl = readGRDECL(grdecl);
    G = processGRDECL(grdecl);
    rock = grdecl2Rock(grdecl, G.cells.indexMap);
    G = computeGeometry(G);
else
    load ~/data/norne.mat
end
%%
close all
set (gcf, 'WindowButtonMotionFcn', @mouseMove);

h = plotCellData(G, rock.poro, 'ButtonDownFcn', @interactiveSelection)
%%
close all;
addFilters(G, [], rock, [])
%%
ijk = gridLogicalIndices(G);
clc;
close all;
data = rock;
plotToolbar2(G, data)
% axis tight off
%%
selection = addFilters(G, {}, data, []);

clf
plotGrid(G, selection)


%%
close all;
extractSubcellsInteractive(G, rock.poro)


%%
mrstModule add spe10
[G, W, rock] = getSPE10setup(25);
rock.poro = max(rock.poro, 1e-4);
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
rS = initState(G, W, 0);
% rS.s(:,1) =
S  = computeMimeticIP(G, rock);
rS = solveIncompFlow(rS, G, S, fluid, 'wells', W);


%% Compute and display time-of-flight and tracer partitioning
% First we compute time-of-flight which is the travel time from an injector
% to a given point in the reservoir, and stationary distribution of tracers
% injected continuously from each injetion well. From the latter, we can
% easily compute the volume flooded by each injector. Reversing the
% velocity field, we can cmopute the reverse time-of-flight (the travel
% time from an arbitrary point to the nearest producer) and the drainage
% volumes of each producer.
D = computeTOFandTracer(rS, G, rock, 'wells', W);
%%
close all
fluidtwo =  initSimpleFluid('mu' , [   1,  10]*centi*poise     , ...
                           'rho', [1014, 859]*kilogram/meter^3, ...
                           'n'  , [   2,   2]);

state = initResSol(G, 0, [.3 .7]);
tmp = rand(G.cells.num, 1);
state.s = [tmp, 1-tmp];


plotTOFArrival(state, fluidtwo, 1, D)
%%



