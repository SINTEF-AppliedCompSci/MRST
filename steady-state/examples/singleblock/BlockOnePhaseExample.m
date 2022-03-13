%% One-Phase Block Upscaling Example
% 
% Example showing permeability upscaling of a single grid block, i.e., the
% entire grid G is upscaled to a single coarse cell. The upscaling can be
% performed using Dirichlet or periodic boundary conditions, or simpler
% averaging methods may be applied.
% 
% Step through the blocks of the script.
% 

%% Add MRST modules

% We rely on the following MRST modules
mrstModule add incomp upscaling steady-state

% In addition, we will use the SPE10 grid as an example model
mrstModule add spe10


%% Construct a model
% We extract a small block from the SPE10 model 1.

% Rock
I = 1:5; J = 30:35; K = 1:5; % Make some selection
rock = getSPE10rock(I, J, K); % Get SPE10 rock

% Grid
cellsize = [20, 10, 2].*ft; % Cell size (ft -> m)
celldim  = [numel(I), numel(J), numel(K)];
physdim  = celldim.*cellsize;
G = cartGrid(celldim, physdim); % Create grid structre
G = computeGeometry(G); % Compute volumes, centroids, etc.


%% Plot the grid

% Create a new figure and set it wide
close(figure(1)); fh = figure(1);
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2.4 1]);

% Permeability on a log-scale
subplot(1,2,1);
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)) );
view(3); axis tight;
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
title('Permeability'); colorbar;

% Porosity
subplot(1,2,2);
plotCellData(G, rock.poro); view(3); axis tight;
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
title('Porosity'); colorbar;


%% Upscale block using Dirichlet BC

% Create GridBlock: to simplify the passing of arguments in different
% upscaling methods, we create a GridBlock instance which holds the grid,
% the rock, and other properties of the grid block that is to be upscaled.
block = GridBlock(G, rock);

% Upscale permeability: upscaled using Dirichlet boundary conditions, and
% no-flow on the normal boundaries. The upscaled data is returned in a
% structure. This structre can grow as more properties are upscaled.
updata = upAbsPermPres(block);

% Upscale porosity of the block by pore volume averaging
updata = upPoro(block, updata);

% Display the upscaled data
updata %#ok<NOPTS>


%% Upscale block using periodic BC

% Create GridBlock with the parameter 'periodic' set to 'true'. The
% GridBlock will then make the grid periodic.
block = GridBlock(G, rock, 'periodic', true);

% Upscale permeability. The function is set up to only return the diagonal
% of the upscaled tensor by default. The previous permeability upscaling is
% replaced.
updata = upAbsPermPres(block, updata) %#ok<NOPTS>

% We can get the full tensor by asking for it.
updata = upAbsPermPres(block, updata, 'fulltensor', true);
updata.perm

% Note that the porosity will not change becuase the grid is periodic.


%% Using the Upscaler class
% 
% Instead of calling the upscaling functions directly, we may use the
% subclasses of the Upscaling class. For one phase upscaling, we call the
% class OnePhaseUpscaler.

% Create an instance of the upscaler 
upscaler = OnePhaseUpscaler(G, rock);

% Perform the upscaling. The data structure returned contains both the
% permeability and the porosity.
updata = upscaler.upscaleBlock(block) %#ok<NASGU,NOPTS>


%% Permeability averaging
% 
% We may also choose to use an averaging method to find uspcaled values of
% the permeability instead of the pressure solver.

% Let us for example compute the arithmetic average.
upscaler.OnePhaseMethod = 'arithmetic';
updata = upscaler.upscaleBlock(block) %#ok<NASGU,NOPTS>

% Another alternative is the combination of harmonic and arithmetic. For
% each different method, observe how the values of the upscaled
% permeability changes.
upscaler.OnePhaseMethod = 'harmonic-arithmetic';
updata = upscaler.upscaleBlock(block) %#ok<NOPTS>

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
