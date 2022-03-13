%% Block Two-Phase Example
% 
% This example demonstrates two-phase upscaling of a single grid block,
% meaning that the entire grid G is upscaled to a single coarse cell. To
% perform a two-phase upscaling, we must first upscale the absolute
% permeability (see <matlab:edit('BlockOnePhaseExample')
% BlockOnePhaseExample>). Subsequently, the relative permeability curves
% are upscaled.
% 
% Step through the code blocks of the script.
% 

%% Add MRST modules

% We rely on the following MRST modules
mrstModule add incomp upscaling ad-props ad-core ad-blackoil steady-state

% In addition, we will use the SPE10 grid as an example model
mrstModule add spe10


%% Open a wide figure

fn = 44; % Just some arbitrary number
close(figure(fn)); fh = figure(fn);
set(fh, 'Name', 'Block Two-Phase Example', 'NumberTitle', 'off');
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
figure(fn); clf; subplot(1,2,1);
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
    'nsat', 30);

% Create an MRST fluid from the property structure
fluid = initADIFluidOW(fprop);


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



%% Create Block and Upscaled Absolute Permeability

% Create GridBlock: to simplify the passing of arguments in different
% upscaling methods, we create a GridBlock instance which holds the grid,
% the rock, and other properties of the grid block that is to be upscaled.
% For two-phase flow, the block also need to access the fluid structure.
block = GridBlock(G, rock, 'fluid', fluid);

% Upscale porosity
updata = upPoro(block);

% Upscale permeability: we need to upscale the (absolute) permeability
% before upscaling relative permeability. Here, we use Dirichlet boundary
% conditions with no-flow on the normal boundaries.
updata = upAbsPermPres(block, updata) %#ok<NOPTS>


%% Relative Permeability upscaling - Viscous-Limit method
% 
% Let's start by upscaling relperm using viscous-limit steady-state
% upscaling. This method assumes the capillary forces are negligable and
% that the fractional flow is constant in the grid block. If nothing is
% specified, then the relative permeability is upscaled in each of the
% three dimensions.

updataVL = upRelPerm(block, updata, 'viscous') %#ok<NOPTS>


%% Plot the upscaled relative permeability
% 
% We plot the upscaled curves with the original curves. Observe that the
% upscaling in x- and y-direction are almost identical, while the
% z-direction upscaling is clearly different. This is not unusual, as
% reservoirs often can have some form of horizontal layering.
% 
% Another point to note, is that for the x- and y-direction, the upscaled
% curves are closer to the original curves of region 1 than region 2. This
% can be explained by region 1 having higher absolute permeability, and
% thus more flow goes through region 1, making this the 'domonant' region.
% 
% Of cource, the above observations will not be valid if the original
% domain is changed.

subplot(1,2,2); cla; hold on;
lh = nan(numel(updataVL.krW)+1,1);
for r=1:2 % Plot original curves in the backgroudn
    lh(1) = plot(fprop.krW{r}(:,1), fprop.krW{r}(:,2), ...
                'Color', [1 1 1].*0.8 );
            plot(1-fprop.krO{r}(:,1), fprop.krO{r}(:,2), ...
                'Color', [1 1 1].*0.8 );
end
colors = lines(numel(updataVL.krW));
for d=1:numel(updataVL.krW) % Plot the upscaled curves on top
    lh(d+1) = plot(updataVL.krW{d}(:,1), updataVL.krW{d}(:,2), ...
                'Color', colors(d,:) );
            plot(updataVL.krO{d}(:,1), updataVL.krO{d}(:,2), ...
                'Color', colors(d,:) );
end
box on; axis([0 1 0 1]); title('Upscaled Relative Permeability');
legend(lh,{'Original','x-dir','y-dir','z-dir'},'Location','North');
xlabel('Water Saturation'); ylabel('Relative Permeability');


%% Relative Permeability upscaling - Capillary-Limit method
% 
% In the other steady-state limit, capillary-limit upscaling, it is assumed
% that the velocity is small enough, such that the capillary forces
% dominate.

% To run the capillary-limit upscaling, we must first upscale the capillary
% pressure curves.
updata = upPcOW(block, updata);

% Then we can upscale the relative permeability
updataCL = upRelPerm(block, updata, 'capillary');


%% Plot the upscaled relative permeabilitis

clf;

for i = 1:2
    subplot(1,2,i); hold on;
    if i==1
        ud = updataVL; % Left: Viscous-limit upscaled curves
        title('Viscous-Limit Upscaling');
    else
        ud = updataCL; % Right: Capillary-limit upscaled curves
        title('Capillary-Limit Upscaling');
    end
    lh = nan(numel(ud.krW)+1,1);
    for r=1:2 % Plot original curves in the backgroudn
        lh(1) = plot(fprop.krW{r}(:,1), fprop.krW{r}(:,2), ...
                    'Color', [1 1 1].*0.8 );
                plot(1-fprop.krO{r}(:,1), fprop.krO{r}(:,2), ...
                    'Color', [1 1 1].*0.8 );
    end
    colors = lines(numel(ud.krW));
    for d=1:numel(ud.krW) % Plot the upscaled curves on top
        lh(d+1) = plot(ud.krW{d}(:,1), ud.krW{d}(:,2), ...
                    'Color', colors(d,:) );
                plot(ud.krO{d}(:,1), ud.krO{d}(:,2), ...
                    'Color', colors(d,:) );
    end
    box on; axis([0 1 0 1]);
    legend(lh,{'Original','x-dir','y-dir','z-dir'},'Location','North');
    xlabel('Water Saturation'); ylabel('Relative Permeability');
end

% Observe how the two different methods produce different upscaled relative
% permeability curves. Especially in the z-direction, the two methods
% differ significantly in this case.



%% Rate-Dependent Steady-State Upscaling

% The viscous- and capillary-limit methods make simplifications which makes
% it faster to find the satuartion distribution at steady state. However,
% in general for intermediate flow-rates, the saturation distribution at
% steady-state must be computed for each value of the average saturation
% and for each flow direction. This flow-based steady-state upscaling is
% *much* more time-consuming, but may be more accurate if none of the
% limits can be assumed.

% Create a block with periodic grid.
blockP = GridBlock(G, rock, 'fluid', fluid, 'periodic', true);

% Run the upscaling. As the steady-state simulations are time-comsuming, we
% supply the verbose option to get updates printed to the console during
% the upscaling.
updataRD = upRelPerm(blockP, updata, 'flow', 'dp', 1*barsa, ...
    'verbose', true);



%% Plot the upscaled relative permeabilitis
% 
% We compare the uscaled curves from the periodic flow-based steady-state
% upscaling with the viscous-limit upscaling.
% 

clf;
for i = 1:2
    subplot(1,2,i); hold on;
    if i==1
        ud = updataVL; % Left: Viscous-limit upscaled curves
        title('Viscous-Limit Upscaling');
    else
        ud = updataRD; % Right: Flow-based upscaling curves
        title('Flow-Based Upscaling');
    end
    lh = nan(numel(ud.krW)+1,1);
    for r=1:2 % Plot original curves in the backgroudn
        lh(1) = plot(fprop.krW{r}(:,1), fprop.krW{r}(:,2), ...
                    'Color', [1 1 1].*0.8 );
                plot(1-fprop.krO{r}(:,1), fprop.krO{r}(:,2), ...
                    'Color', [1 1 1].*0.8 );
    end
    colors = lines(numel(ud.krW));
    for d=1:numel(ud.krW) % Plot the upscaled curves on top
        lh(d+1) = plot(ud.krW{d}(:,1), ud.krW{d}(:,2), ...
                    'Color', colors(d,:) );
                plot(ud.krO{d}(:,1), ud.krO{d}(:,2), ...
                    'Color', colors(d,:) );
    end
    box on; axis([0 1 0 1]);
    legend(lh,{'Original','x-dir','y-dir','z-dir'},'Location','North');
    xlabel('Water Saturation'); ylabel('Relative Permeability');
end

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
