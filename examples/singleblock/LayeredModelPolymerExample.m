%% Layered Model Polymer Upscaling Example
% 
% This example demonstrates polymer upscaling of a simple layered model.
% The upscaling of polymer is described in [1], and this example is similar
% to an example in that paper. There exists an analytical solution to the
% upscaling, and so this example serves as a verification.
% 
% Step through the code blocks of the script.
% 
% References:
% [1] Hilden, Lie and Xavier Raynaud - "Steady State Upscaling of Polymer
%     Flooding", ECMOR XIV - 14th European conference on the mathematics of
%     oil recovery, 2014.
% 


%% Add MRST modules

% Add this module
mrstModule add steady-state

% We rely on the following MRST modules
mrstModule add incomp upscaling ad-props ad-core ad-blackoil ad-fi

%% Open a new figure and make it wide

fh = figure;
set(fh, 'Name', 'Layered Model Example', 'NumberTitle', 'off');
op = get(fh, 'OuterPosition');
set(fh, 'OuterPosition', op.*[1 1 2.4 1]);

%% Create grid

% Grid structure
G = cartGrid([3 3 3]);
G = computeGeometry(G);
ijk = gridLogicalIndices(G);

% Set middel layer is rock #2, while top and bottom layers are rock #1
regnum = ones(G.cells.num,1);
regnum(ijk{3}==2) = 2; 

% Rock
K = [100 0.1].*milli*darcy;
rock.perm = K(1).*ones(G.cells.num,1);
rock.perm(regnum==2) = K(2);
rock.poro = 0.1.*ones(G.cells.num,1);

% Fluid with different relperms in the two regions
fprop = getExampleFluidProps(rock, 'satnum', regnum, ...
    'swir', 0, 'sor', 0, 'krWmax', 1, 'krn', [2.5 1.5], ...
    'nsat', 50, 'polymer', true);
fprop.pcOW = [];
fprop.rhoO = 1000;
fprop.rhoW = 1000;
fprop.muO = 1*centi*poise;
fprop.muW = 1*centi*poise;

% Polymer properties
fprop.ads = {[ [ 0; 0.3; 1.5; 2.6;  4].*kilogram/meter^3, ...
               [ 0;  15;  27;  30; 30].*(milli*gram)/(kilo*gram) ], ...
             [ [ 0; 0.3; 1.5; 2.7;  4].*kilogram/meter^3, ...
               [ 0;  10;  20;  25; 25].*(milli*gram)/(kilo*gram) ] };
fprop.adsMax = [30; 25].*(milli*gram)/(kilo*gram);
fprop.muWMult = [ [ 0;  1;  3;  4].*kilogram/meter^3, ...
                  [ 1;  3; 37; 40] ];
fprop.rrf = [3.2; 1.4];

Rkorg = cell(1,2);
for r=1:2
    Rkorg{r} = [fprop.ads{r}(:,1), ...
        1 + (fprop.rrf(r)-1).*( fprop.ads{r}(:,2)./fprop.adsMax(r) ) ];
end

fluid = initADIFluidOWPolymer(fprop);


%% Plot the grid

% Plot the regions
clf; subplot(1,2,1);
plotCellData(G, regnum);
view(3); axis tight; cbh = colorbar; title('Rock types');
set(cbh, 'YTick', [1 2]);

% Plot the permeability
subplot(1,2,2);
plotCellData(G, log10(convertTo(rock.perm(:,1), milli*darcy)));
view(3); axis tight; colorbar; title('Permeability (mD)');


%% Plot relative permeabilities

% Plot the regions
clf; subplot(1,2,1);
plotCellData(G, regnum);
view(3); axis tight; cbh = colorbar; title('Rock types');
set(cbh, 'YTick', [1 2]);

% Plot the relative permeability
subplot(1,2,2); cla; hold on;
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


%% Plot polymer properties

% Adsorption
clf; subplot(1,2,1); hold on;
colors = lines(2); lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.ads{r}(:,1), fprop.ads{r}(:,2).*(10^6), ...
                'Color', colors(r,:) );
end
box on; title('Adsorption'); axis([0 4 0 35]);
legend(lh,{'Region 1','Region 2'},'Location','NorthWest');
xlabel('Polymer Concentration (kg/m^3)'); ylabel('Adsorption (mg/kg)');

% Reduction Factor
subplot(1,2,2); hold on;
lh = nan(2,1);
for r=1:2
    lh(r) = plot(fprop.ads{r}(:,2).*(10^6), Rkorg{r}(:,2));
end
legend(lh,{'Region 1','Region 2'},'Location','NorthWest');
box on; title('Reduction Factor, Rk'); axis([0 30 1 3.5]);
xlabel('Adsorption (mg/kg)'); ylabel('Reduction Factor, Rk');


%% Upscale absperm

% Create grid block
block = GridBlock(G, rock, 'fluid', fluid, 'periodic', true);

% Upscale porosity
updata = upPoro(block);

% Upscale absolute permeability
updata = upAbsPerm(block, updata);

% Get analytical solution
updataA = LayeredModelExactUpscaling(K, fprop.krW, fprop.krO, Rkorg);


%% Absolute permeability
% 
% To verify that the absolute permeability is correctly upscaled, we
% compare the numerically computed values with the analytical solution.
% 

KA = updataA.perm./(milli*darcy);
KN = updata.perm./(milli*darcy);
fprintf('\nUpscaled absolute permeability\n');
fprintf('Analytical: [%1.2f, %1.2f, %1.2f]\n', KA(1), KA(2), KA(3));
fprintf('Numerical:  [%1.2f, %1.2f, %1.2f]\n', KN(1), KN(2), KN(3));
fprintf('Rel.diff.:  %1.2e\n\n', norm( (KA-KN)./KA ) );


%% Upscale relative permeability and polymer reduction factor Rk
% 
% We upscale the relative permeability using flow-based steady-state
% upscaling with periodic boundary conditions (periodicity is determined by
% the creation of the block structure above).
% 
% The flow-based upscaling is rather time-consuming, as a two-phase flow
% simulation is run for each upscaled saturation value and each dimension.
% Requesting 15 upscaled values, this amounts to 13*3=39 flow simulations.
% We do not need to run the endpoints, as one of the phases is then
% immobile and the solution is known a priori.
% 
% The polymer reduction factor is upscaled simultaneously as the relative
% permeability. At steady state, the polymer concentration is constant in
% the domain (see [1]). Therefore, the saturation distribution do not need
% to be recomputed for the polymer upscaling, but we reuse the distribution
% from the relative permeability upscaling.
% 

% Upscale relative permeability AND reduction factor Rk
updata = upRk(block, updata, 'flow', 'nsat', 15, 'npoly', 10, ...
                'verbose', true);

% Upscale polymer adsorption isoterm
updata = upPolyAds(block, updata);

% Note: If we only wanted to upscale the relative permeability, we wound
% instead make the call
%    updata = upRelPerm(block, updata, ...);


%% Plot upscaled relative permeability
% 
% The upscaling of the relative permeabilities are also computed
% analytically, and we compare the analytical solution with the numerically
% flow-based upscaled curves. This plot is equivalent to Figure 2 in [1].
% 

clf;
lw = 1;
colors = lines(numel(updata.krW));
lh = nan(4,1);
dims = [1 3];

% Plot both x- and z-direction. y is equal to x.
for i=1:2
    d = dims(i);
    subplot(1,2,i); hold on;
    
    % Plot original curves first
    for r=1:2
        lh(r) = plot(fprop.krW{r}(:,1), fprop.krW{r}(:,2), ...
                    'Color', colors(r,:), 'LineWidth', lw );
                plot(1-fprop.krO{r}(:,1), fprop.krO{r}(:,2), ...
                    'Color', colors(r,:), 'LineWidth', lw );
    end
    
    % Plot analytical upscaling and numerical upscaling
    lh(3) = plot(updataA.krW{d}(:,1), updataA.krW{d}(:,2), 'k--', ...
        'LineWidth', lw );
    plot(updataA.krO{d}(:,1), updataA.krO{d}(:,2), 'k--', ...
        'LineWidth', lw );
    lh(4) = plot(updata.krW{d}(:,1), updata.krW{d}(:,2), 'ko', ...
        'LineWidth', lw );
    plot(updata.krO{d}(:,1), updata.krO{d}(:,2), 'ko', ...
        'LineWidth', lw );
    
    box on; axis([0 1 0 1]);
    if (i==1),title('x/y-direction');else title('z-direction');end
    legend(lh,{'Rock 1','Rock 2','Analytical','Numerical'}, ...
        'Location','North');
    xlabel('Water Saturation'); ylabel('Relative Permeability');
end

%% Plot upscaled reduction factor Rk
% 
% 
% 

clf;
lw = 1;
colors = lines(numel(updata.krW));
lh = nan(4,1);
dims = [1 3];

% Plot both x- and z-direction. y is equal to x.
for i=1:2
    d = dims(i);
    subplot(1,2,i); hold on;
    
    % Plot original curves first
    for r=1:2
        lh(r) = plot(fprop.ads{r}(:,2).*(10^6), Rkorg{r}(:,2), ...
                    'Color', colors(r,:), 'LineWidth', lw );
    end
    
    % Plot numerical upscaling
    c   = updata.Rk.c;
    ads = interp1(updata.ads(:,1), updata.ads(:,2), c);
    sW  = updata.Rk.s{d};
    for is=1:numel(sW)
        fprintf('sW=%1.2f\n', sW(is));
        plot(ads.*(10^6), updata.Rk.val{d}(is,:), 'k' );
    end
    
%     lh(3) = plot(updataA.krW{d}(:,1), updataA.krW{d}(:,2), 'k--', ...
%         'LineWidth', lw );
%     plot(updataA.krO{d}(:,1), updataA.krO{d}(:,2), 'k--', ...
%         'LineWidth', lw );
%     lh(4) = plot(updata.krW{d}(:,1), updata.krW{d}(:,2), 'ko', ...
%         'LineWidth', lw );
%     plot(updata.krO{d}(:,1), updata.krO{d}(:,2), 'ko', ...
%         'LineWidth', lw );
    
    if (i==1),title('x/y-direction');else title('z-direction');end
%     legend(lh,{'Rock 1','Rock 2','Analytical','Numerical'}, ...
%         'Location','North');
    box on; title('Reduction Factor, Rk'); axis([0 30 1 3.5]);
    xlabel('Adsorption (mg/kg)'); ylabel('Reduction Factor, Rk');
end


