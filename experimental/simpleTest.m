%% Simple test of the upscaling

clear variables

G = cartGrid([4,4,2]);
G = computeGeometry(G);

mrstVerbose on

rock.perm = rand(G.cells.num,1).*(500*milli*darcy);
rock.poro = rand(G.cells.num,1).*0.3 + 0.1;

fluid = initSimpleFluid(...
    'mu' , [   1,  10]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   2,   2]);

upscaler = OnePhaseUpscaler(G, rock, fluid);



cdims = [2 2 1];

% Partition grid
p = partitionUI(G, cdims);
p = compressPartition(p);

[data, report] = upscaler.upscale(p);

