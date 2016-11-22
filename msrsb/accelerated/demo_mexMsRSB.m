mrstModule add coarsegrid msrsb libgeometry incomp
% Some different test cases
testcase = 'simple';

switch lower(testcase)
    case 'simple'
        G = cartGrid([15, 10]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [3, 2]);
    case 'medium'
        G = cartGrid([80, 40]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [8, 4]);
    case 'big'
        G = cartGrid([50, 50, 20]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [5, 5, 2]);
    case 'huge'
        G = cartGrid([200, 200, 100]);
        rock = makeRock(G, 1, 1);
        p = partitionUI(G, [10, 10, 10]);
        G = mcomputeGeometry(G);
    case 'spe10'
        mrstModule add spe10
        layers = 1:85;
        [G, ~, rock] = getSPE10setup(layers);
        p = partitionUI(G, [6, 11, ceil(G.cartDims(3)./5)]);
    otherwise
        error();
end

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);

CG = setupMexInteractionMapping(CG);
%% Generate basis functions via regular mex interface
% Typical usage. We send in the coarsegrid and the mex layer passes values
% onto the C++ code.
I = cppMultiscaleBasis(CG, A, 'verbose', true, 'omega', 2/3, 'maxiter', 250, 'tolerance', 0.001);

%% Write a testcase to flatfiles on disk
% We can also write basis functions to disk as a series of text files. This
% is not efficient, but it is useful for setting up and running testcases
% without Matlab.
fn = fullfile(mrstDataDirectory(), 'tmp', ['basis_', lower(testcase)]);
if ~exist(fn, 'dir')
    cppMultiscaleBasis(CG, A, 'doSolve', false, 'writePath', fn);
end

disp('Stored files:')
ls(fullfile(fn, 'input'))
%% Pass the filename to the mex function
% The function will read the files and produce output in a folder named
% output. Once it is done, it will be read in and returned as the same type
% of sparse matrix.
I2 = cppMultiscaleBasisFromFile(fn);

%% Plot the basis functions
mrstModule add mrst-gui
close all;

figure;
plotToolbar(G, I);
colorbar
title('Basis functions made from mex gateway')

figure;
plotToolbar(G, I2);
colorbar
title('Basis functions made from txt files')
