%% Block Polymer Example
% 
% This example demonstrates polymer upscaling of a single grid block,
% meaning that the entire grid G is upscaled to a single coarse cell.
% 
% The upscaling of polymer is described in [1]. First, we upscale the
% absolute permeability and the relative permeavbilities. Then, the polymer
% reduction factor Rk is upscaled. The upscaling methodology resebles the
% upscaling of relative permeability.
% 
% Step through the code blocks of the script.
% 
% References:
% [1] Hilden, Lie and Xavier Raynaud - "Steady State Upscaling of Polymer
%     Flooding", ECMOR XIV - 14th European conference on the mathematics of
%     oil recovery, 2014.
% 

%% Add MRST modules

% We rely on the following MRST modules
mrstModule add incomp upscaling ad-props ad-core ad-blackoil steady-state

% In addition, we will use the SPE10 grid as an example model
mrstModule add spe10


%% Open a new figure and make it wide

fh = figure;
set(fh, 'Name', 'Block Polymer Example', 'NumberTitle', 'off');
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2.4 1]);


%% Construct a model
% We extract a small block from the SPE10 model 1.

% Rock
I = 1:5; J = 30:35; K = 1:5; % Make some selection
rock = getSPE10rock(I, J, K); % Get SPE10 rock
rock.poro(rock.poro==0) = min(rock.poro(rock.poro>0)); % Remove zero poro

% Grid
cellsize = [20, 10, 2].*ft; % Cell size (ft -> m)
celldim  = [numel(I), numel(J), numel(K)];
physdim  = celldim.*cellsize;
G = cartGrid(celldim, physdim); % Create grid structre
G = computeGeometry(G); % Compute volumes, centroids, etc.


%% Plot the grid

% Plot the permeability
clf; subplot(1,2,1);
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)));
view(3); axis tight; colorbar; title('Permeability');

% Plot the regions
subplot(1,2,2);
plotCellData(G, rock.poro); view(3); axis tight;
colorbar; title('Porosity');


%% Create two regions

% Create two regions which will have different relative permeability data.
% This will make the relative permeability upscaling more interesting. We
% divide the cells in two equal sizes, based on the value of the
% permeability.
kl = log10(rock.perm(:,1));
regnum = ones(G.cells.num, 1);
regnum(kl<median(kl)) = 2;



%% Plot the permeability and regions

% Plot the permeability
clf; subplot(1,2,1);
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)));
view(3); axis tight; colorbar; title('Permeability');

% Plot the regions
subplot(1,2,2);
plotCellData(G, regnum); view(3); axis tight;
colorbar; title('Rock regions');


%% Create a fluid
% We apply some helper functions to help us create an example fluid

% Get a property structure
fprop = getExampleFluidProps(rock, 'satnum', regnum, ...
    'swir', [0.1 0.25], 'sor',[0.14 0.18], 'krWmax', [0.8 0.6], ...
    'nsat', 30, 'polymer', true);

% Create an MRST fluid from the property structure
fluid = initADIFluidOWPolymer(fprop);


%% Plot fluid properties

clf;

% Relative permeability
subplot(1,2,1); hold on;
colors = lines(2); lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.krW{r}(:,1), fprop.krW{r}(:,2), ...
                'Color', colors(r,:) );
            plot(1-fprop.krO{r}(:,1), fprop.krO{r}(:,2), ...
                'Color', colors(r,:) );
end
box on; axis([0 1 0 1]); title('Relative Permeability');
legend(lh,{'Region 1','Region 2'},'Location','North');
xlabel('Water Saturation'); ylabel('Relative Permeability');

% Capillary pressure
subplot(1,2,2); hold on;
lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.pcOW{r}(:,1), fprop.pcOW{r}(:,2)./barsa, ...
                'Color', colors(r,:) );
end
box on; axis([0 1 -1 1]); title('Capillary Pressure');
legend(lh,{'Region 1','Region 2'},'Location','North');
xlabel('Water Saturation'); ylabel('Capillary Pressure (bar)');


%% Plot polymer properties

clf;

% Relative permeability
subplot(1,2,1); hold on;
colors = lines(2); lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.ads{r}(:,1), fprop.ads{r}(:,2).*(10^6), ...
                'Color', colors(r,:) );
end
box on; title('Adsorption'); axis([0 4 0 80]);
legend(lh,{'Region 1','Region 2'},'Location','NorthWest');
xlabel('Polymer Concentration (kg/m^3)'); ylabel('Adsorption (mg/kg)');

% Capillary pressure
subplot(1,2,2); hold on;
lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.muWMult{r}(:,1), fprop.muWMult{r}(:,2), ...
                'Color', colors(r,:) );
end
box on; title('Viscosity Multiplier'); axis([0 4 0 45]);
legend(lh,{'Region 1','Region 2'},'Location','NorthWest');
xlabel('Polymer Concentration (kg/m^3)'); ylabel('Viscosity Multiplier');


%% Create Block and Upscale Absolute and Relative Permeability
% 
% See the exampels BlockOnePhaseExample and BlockTwoPhaseExample for more
% on one- and two-phase upscaling.
% 

% Create grid block
block = GridBlock(G, rock, 'fluid', fluid);

% Upscale porosity
updata = upPoro(block);

% Upscale absolute permeability
updata = upAbsPermPres(block, updata);

% Upscale relative permeability (using viscous-limit method)
%updata = upRelPerm(block, updata, 'viscous');


%% Upscale Permeability

updata = upRelPerm(block, updata, 'viscous', 'dims', 1);

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
