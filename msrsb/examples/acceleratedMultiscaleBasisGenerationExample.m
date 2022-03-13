%% MsRSB basis function with mex acceleration
% The MsRSB basis functions can be efficiently computed using a compiled
% function written in C++. If a C++ compiler is available, this is
% generally an order of magnitude faster than using the Matlab version.
%
% In this example, we demonstrate the syntax for generating basis functions
% either directly from Matlab or via a flat file format.
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
        error('Unknown test case %s', testcase);
end

if ~isfield(G.cells, 'centroids')
    G = computeGeometry(G);
end

T = computeTrans(G, rock);
A = getIncomp1PhMatrix(G, T);

CG = generateCoarseGrid(G, p);
CG = coarsenGeometry(CG);
CG = storeInteractionRegionCart(CG);
%% Add extra data to coarse grid
% The function setupMexInteractionMapping adds additional mappings to the
% coarse grid which are required for the compiled interface. This can be
% done once for a specific coarse-fine grid pair and reused throughout a
% simulation.
CG = setupMexInteractionMapping(CG);
%% Generate basis functions via regular mex interface
% Typical usage. We send in the coarsegrid and the mex layer passes values
% onto the C++ code.
I = cppMultiscaleBasis(CG, A, 'verbose', true, 'omega', 2/3,...
                              'maxiter', 250, 'tolerance', 0.001);
%% Generate basis functions via general interface
% This is also supported through the general interface for basis functions.
basis = getMultiscaleBasis(CG, A, 'useMex', true, 'type', 'msrsb');
%% Write a testcase to flatfiles on disk
% We can also write basis functions to disk as a series of text files. This
% is not efficient, but it is useful for setting up and running testcases
% without running a full Matlab session, or for generating basis functions
% for some other purpose.
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

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
