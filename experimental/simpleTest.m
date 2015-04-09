%% Simple test of the upscaling

mrstModule add coarsegrid upscaling
mrstModule add deckformat ad-fi

clear variables

fn = '/home/sindre/git/mrst-my-modules/projects/fieldModel/models/temp/fine/MODEL.DATA';



deck  = readEclipseDeck(fn);
deck  = convertDeckUnits(deck);
G     = initEclipseGrid(deck);
G     = computeGeometry(G);
rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);
fluid = initDeckADIFluid(deck);

% Add fields for two phase upscaling
fluid = addPcOWInvADIFluid(fluid, deck);
fluid = addFracFlowInvADIFluid(fluid, deck);

% G = cartGrid([4,4,2]);
% G = computeGeometry(G);

% mrstVerbose on

% rock.perm = rand(G.cells.num,1).*(500*milli*darcy);
% rock.poro = rand(G.cells.num,1).*0.3 + 0.1;

% fluid = initSimpleFluid(...
%     'mu' , [   1,  10]*centi*poise     , ...
%     'rho', [1014, 859]*kilogram/meter^3, ...
%     'n'  , [   2,   2]);

% cdims = [2 2 1];

cdims = [2 2 1];

% Partition grid
p = partitionUI(G, cdims);
p = compressPartition(p);

%%

figure;
plotCellData(G, p);
view(3); axis tight;

%%

mrstVerbose off
t = tic;
upscaler = OnePhaseUpscaler(G, rock);
upscaler.partition = p;
%upscaler.blocks = 2;
upscaler.dims = 1:3;
dataOP = upscaler.upscale();
toc(t)


%%



mrstVerbose off
t = tic;
upscaler = TwoPhaseUpscaler(G, rock, fluid);
upscaler.partition = p;
%upscaler.blocks = 2;
upscaler.dims = 1;
upscaler.nvalues = 19;
upscaler.method = 'capillary';
upscaler.verbose = true;
dataTP = upscaler.upscale();
toc(t)




%% Old method

% Create coarse grid
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG.cartDims = cdims;
CG.cells.indexMap = (1:CG.cells.num)';

mrstVerbose off
t = tic;
upscaled = exampleUpscale(G, rock, fluid, 'CG', CG, ...
    'method', 'capillary', 'save', false, 'periodic', false, ...
    'dims', 1);
toc(t)

% 
% % KU = upAbsPermADI(G, CG, rock, 'periodic',false);
% 

%% Check that the upscaling is the same

for i=1:max(p)
    
    norm( data(i).K(:) - KU(i,:)' )
    
end


