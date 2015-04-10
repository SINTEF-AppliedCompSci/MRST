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

% figure;
% plotCellData(G, p);
% view(3); axis tight;

%% Settings

dims = 1;
%method = 'capillary';
method = 'viscous';

%%



mrstVerbose off
t = tic;
upscaler = TwoPhaseUpscaler(G, rock, fluid);
upscaler.partition = p;
%upscaler.blocks = 2;
upscaler.dims     = dims;
upscaler.nrelperm = 19;
upscaler.npcow    = 50;
upscaler.method   = method;
upscaler.verbose  = true;
dataTP = upscaler.upscale();
toc(t)

%%

figure; hold on;
arrayfun(@(x) plot(x.krW(:,1), x.krW(:,2)), dataTP);
arrayfun(@(x) plot(x.krO(:,1), x.krO(:,2)), dataTP);


%%



mrstVerbose off
t = tic;
upscaler = TwoPhaseStepwiseUpscaler(G, rock, fluid);
upscaler.partition = p;
upscaler.nrelperm  = 19;
upscaler.npcow     = 200;
upscaler.method1   = 'capillary';
upscaler.dim1      = 3;
upscaler.method2   = 'viscous';
upscaler.dim2      = 1;
upscaler.verbose   = true;
dataStep = upscaler.upscale();
toc(t)

%%

figure; hold on;
arrayfun(@(x) plot(x.krW(:,1), x.krW(:,2)), dataStep);
arrayfun(@(x) plot(x.krO(:,1), x.krO(:,2)), dataStep);


%% Old method

% Create coarse grid
CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG.cartDims = cdims;
CG.cells.indexMap = (1:CG.cells.num)';

mrstVerbose off
t = tic;
upscaled = exampleUpscale(G, rock, fluid, 'CG', CG, ...
    'method', method, 'save', false, 'periodic', false, ...
    'dims', dims);
toc(t)

% 
% % KU = upAbsPermADI(G, CG, rock, 'periodic',false);
% 

%% Check that the upscaling is the same

maxNorm = 0;
for i=1:max(p)
    n = norm( dataTP(i).K(:) - upscaled.perm(i,:)' );
    maxNorm = max(maxNorm, n );
end
fprintf('Max norm K:   %1.2e\n', maxNorm);

maxNorm = 0;
for i=1:max(p)
    for j=1:numel(dims)
        n = norm(dataTP(i).krW{j}(:,2)-upscaled.krRaw{i}.krW(:,j));
        maxNorm = max(maxNorm, n);
    end
end
fprintf('Max norm krW: %1.2e\n', maxNorm);

maxNorm = 0;
for i=1:max(p)
    for j=1:numel(dims)
        n = norm(dataTP(i).krO{j}(:,2)-upscaled.krRaw{i}.krO(:,j));
        maxNorm = max(maxNorm, n);
    end
end
fprintf('Max norm krO: %1.2e\n', maxNorm);

maxNorm = 0;
for i=1:max(p)
    n = norm(dataTP(i).pcOW(:,2) - upscaled.pcOW{i}(:,2));
    maxNorm = max(maxNorm, n);
end
fprintf('Max norm pc:  %1.2e\n', maxNorm);



