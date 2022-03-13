%% Layered Model Polymer Upscaling Example
% 
% This example demonstrates polymer upscaling of a simple layered model.
% The upscaling of polymer is described in [1], and this example is
% identical the first example in therein. There exists an analytical
% solution to the upscaling, and so this example serves as a verification.
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
mrstModule add incomp upscaling ad-props ad-core ad-blackoil

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
clear rock
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
legend(lh,{'Rock 1','Rock 2'},'Location','North');
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
legend(lh,{'Rock 1','Rock 2'},'Location','NorthWest');
box on; title('Reduction Factor, Rk'); axis([0 30 1 3.5]);
xlabel('Adsorption (mg/kg)'); ylabel('Reduction Factor, Rk');


%% Upscale absperm

% Create grid block
block = GridBlock(G, rock, 'fluid', fluid, 'periodic', true);

% Upscale porosity
updata = upPoro(block);

% Upscale absolute permeability
updata = upAbsPerm(block, updata);


%% Compute analytical solution
% 
% This example has an analytical solution for the absolute permeability,
% relative permeability and the reduction factor, as described in [1]. This
% allows us to compare the numerical flow-based solution with the
% analytical solution as a verification.
% 

updataA = LayeredExactUpscaling(K, fprop.krW, fprop.krO, Rkorg);


%% Compare absolute permeability upscaling
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


%% Upscale relative permeability
% 
% We upscale the relative permeability using flow-based steady-state
% upscaling with periodic boundary conditions (periodicity is determined by
% the creation of the block structure above).
% 
% The flow-based upscaling is rather time-consuming, as a two-phase flow
% simulation is run for each upscaled saturation value and each dimension.
% As the model is symmetric in x-y, we only upscale in the x- and
% z-direction to save computation. But still, e.g. requesting 15 upscaled
% values, the number of flow simulations amounts to 13*2=26 . We do not
% need to run the endpoints, as one of the phases is then immobile and the
% solution is known a priori.
% 

% Upscale relative permeability
updata = upRelPerm(block, updata, 'flow', 'nsat', 15, 'dims', [1 3], ...
    'verbose', true);



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
    
    % Plot analytical upscaling
    lh(3) = plot(updataA.krW{d}(:,1), updataA.krW{d}(:,2), 'k--', ...
        'LineWidth', lw );
    plot(updataA.krO{d}(:,1), updataA.krO{d}(:,2), 'k--', ...
        'LineWidth', lw );
    
    % Plot numerical upscaling. Note that the second dimension is z, so we
    % use index i instead of d.
    lh(4) = plot(updata.krW{i}(:,1), updata.krW{i}(:,2), 'ko', ...
        'LineWidth', lw );
    plot(updata.krO{i}(:,1), updata.krO{i}(:,2), 'ko', ...
        'LineWidth', lw );
    
    box on; axis([0 1 0 1]);
    if (i==1),title('x/y-direction');else title('z-direction');end
    legend(lh,{'Rock 1','Rock 2','Analytical','Numerical'}, ...
        'Location','North');
    xlabel('Water Saturation'); ylabel('Relative Permeability');
end


%% Polymer upscaling
% 
% The polymer upscaling is performed as described in [1]. The adsorption
% isoterm is obtained by a simple average, while the reduction factor Rk is
% upscaled in a similar fashion as the relative permeability. However, this
% upscaling may depend on both water saturation and polymer concentration,
% and the upscaling cost becomes even larger.
% 

% Upscale polymer adsorption isoterm
updata = upPolyAds(block, updata);

% Upscale polymer reduction factor Rk
updata = upPolyRk(block, updata, 'flow', 'nsat', 5, 'npoly', 5, ...
    'dims', [1 3], 'verbose', true);



%% Plot upscaled reduction factor Rk
% 
% Finally, we compare the upscaled values of Rk. The plot generated below
% is equivalent to Figure 3 in [1]. Notice the saturation dependency of the
% upscaled Rk in the z-direction.
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
    
    % Plot analytical upscaling
    c   = updataA.Rk.c;
    ads = interp1(updata.ads(:,1), updata.ads(:,2), c);
    sW  = updata.Rk.s{i}; % Numerical
    for is=1:numel(sW)
        s  = sW(is);
        Rk = nan(1,numel(c));
        for ic=1:numel(c)
            Rk(ic) = interp1(updataA.Rk.s{d}, updataA.Rk.val{d}(:,ic), s);
        end
        lh(3) = plot(ads.*(10^6), Rk, 'k', 'LineWidth', lw );
    end
    
    % Plot numerical upscaling. Note that the second dimension is z, so we
    % use index i instead of d.
    c   = updata.Rk.c;
    ads = interp1(updata.ads(:,1), updata.ads(:,2), c);
    sW  = updata.Rk.s{i};
    for is=1:numel(sW)
        lh(4) = plot(ads.*(10^6), updata.Rk.val{i}(is,:), 'ko', ...
            'LineWidth', 1.5 );
    end
    
    if (i==1),title('x/y-direction');else title('z-direction');end
    legend(lh, {'Rock 1','Rock 2','Analytical','Numerical'}, ...
        'Location','NorthWest');
    box on; axis([0 30 1 3.5]);
    xlabel('Adsorption (mg/kg)'); ylabel('Reduction Factor, Rk');
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
