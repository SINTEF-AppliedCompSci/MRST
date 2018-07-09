%% Two-Phase Problem with a Quarter Five-Spot Well Pattern
% In this second introductory example to the HFM module, we show the impact
% of fractures on fluid migration using the hierarchical/embedded fracture
% model. To this end, we consider a two-phase example with three
% intersecting fractures in the center of the model. Oil is recovered by a
% production well in the NE corner, which is supported by a water-injector
% in the SW corner. 

% Load necessary modules, etc 
mrstModule add hfm;             % hybrid fracture module
mrstModule add mrst-gui;        % plotting routines
mrstModule add incomp;          % Incompressible fluid models
checkLineSegmentIntersect;      % ensure lineSegmentIntersect.m is on path

%% Grid and fracture lines
% Construct a Cartesian grid comprising 50-by-20 cells, where each cell has
% dimension 10-by-10 m^2. Define 3 fracture lines by supplying their end
% points in the form [x1 y1 x2 y2].

celldim = [50 20];
physdim = [500 200];
G = cartGrid(celldim, physdim);
G = computeGeometry(G);

fl = [80,  100, 270, 180;...
      130, 160, 340,  40;...
      260,  40, 420, 150] ; % fractures lines in [x1 y1 x2 y2] format.

%% Process fracture lines
% Using the input fracture lines, we identify independent fracture networks
% comprising of connected lines. In this example, there is only 1 fracture
% network consisting of 3 fracture lines. We also identify the fine-cells
% in the matrix containing these fractures. Fracture aperture is set to
% 0.04 meters. The matrix grid and fracture lines are plotted.

dispif(mrstVerbose, 'Processing user input...\n\n');
[G,fracture] = processFracture2D(G, fl, 'verbose', mrstVerbose);
fracture.aperture = 1/25; % Fracture aperture
figure;
plotFractureLines(G,fracture);
axis equal tight; 
box on

%% Compute CI and construct fracture grid
% For each matrix block containing a fracture, we compute a fracture-matrix
% conductivity index (CI) in accordance with the hierarchical fracture
% model (a.k.a. embedded discrete fracture model). Following that, we
% compute a fracture grid where the fracture cell size is defined to be
% 10 m. Next, the fracture grid is plotted on top of the matrix grid.

dispif(mrstVerbose, 'Computing CI and constructing fracture grid...\n\n');
G = CIcalculator2D(G,fracture);
min_size = 5; cell_size = 10; % minimum and average cell size.
[G,F,fracture] = gridFracture2D(G,fracture,'min_size',min_size,'cell_size',cell_size);
clf; plotFractureNodes2D(G,F,fracture); 
axis equal tight; box on

%% Set rock properties in fracture and matrix
% For this example, we will generate the porosity as a Gaussian field. To
% get a crude approximation to the permeability-porosity relationship, we
% assume that our medium is made up of uniform spherical grains of diameter
% dp = 10 m, for which the specic surface area is Av = 6 = dp. With these
% assumptions, using the Carman Kozeny relation, we can then calculate the
% isotropic permeability (K). The rock properties are then plotted.

dispif(mrstVerbose, 'Initializing rock and fluid properties...\n\n');
p = gaussianField(celldim, [0.2 0.4], [11 3], 2.5);
K = p.^3.*(1e-5)^2./(0.81*72*(1-p).^2);
G.rock.poro = p(:);
G.rock.perm = K(:);
K_frac = 10000; % Darcy
G = makeRockFrac(G, K_frac, 'porosity', 0.8);
clf; plotToolbar(G, G.rock); axis equal tight;
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);

%% Define fluid properties
% Define a two-phase fluid model without capillarity. The fluid model has
% the values:
%
% * densities: [rho_w, rho_o] = [1000 700] kg/m^3
% * viscosities: [mu_w, mu_o] = [1 1] cP.
% * corey-coefficient: [2, 2] = [2 2].

fluid = initSimpleFluid('mu' , [   1,  1] .* centi*poise     , ...
    'rho', [1000, 700] .* kilogram/meter^3, ...
    'n'  , [   2,   2]);

%% Define fracture connections as NNC and compute the transmissibilities
% In this section, we use the function defineNNCandTrans to combine the
% fracture and matrix grid structures into a single grid structure. In
% addition to that, we assign a 'non-neighbouring connection (NNC)' status
% to every fracture-matrix connection. For example, if fracture cell 'i' is
% located in the matrix cell 'j', then the connection i-j is defined as an
% NNC. To compute the flux between these elements, we compute a
% transmissibility for each NNC using the CI's computed earlier. Vector T
% contains the transmissibility for each face in the combined grid and each
% NNC.

[G,T] = defineNNCandTrans(G,F,fracture);

%% Add wells
% We will include two wells, one rate-controlled injector well and the
% producer controlled by bottom-hole pressure. The injector and producer
% are located at the bottom left and top right corner of the grid,
% respectively. Wells are described using a Peaceman model, giving an extra
% set of equations that need to be assembled, see simpleWellExample.m for
% more details on the Peaceman well model.

inj = 1;
prod = celldim(1)*celldim(2);
wellRadius = 0.1;

W = addWell([], G.Matrix, G.Matrix.rock, inj,'InnerProduct', 'ip_tpf','Type', ...
    'rate', 'Val', 10, 'Radius', wellRadius, 'Comp_i', [1, 0]);
W = addWell(W, G.Matrix, G.Matrix.rock, prod, 'InnerProduct', 'ip_tpf', 'Type', ...
    'bhp' , 'Val', 50*barsa, 'Radius', wellRadius, 'Comp_i', [0, 1]);

%% Initialize state variables
% Here, we initialize the solution structure using the combined grid and
% the wells defined above.

dispif(mrstVerbose, 'Computing coefficient matrix...\n\n');
state  = initResSol (G, 0);
state.wellSol = initWellSol(W, 0);

%% Compute initial pressure

dispif(mrstVerbose, '\nSolving fine-scale system...\n\n');
state = incompTPFA(state, G, T, fluid, 'Wells', W, 'MatrixOutput', true, 'use_trans',true);

%% Plot initial pressure
clf
hp = plotCellData(G, state.pressure,'EdgeColor','none');
xw = G.cells.centroids(vertcat(W.cells),:);
hold on
plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
hold off
colormap jet
view(0, 90); colorbar; axis equal tight;
title('Initial pressure');

%% Incompressible Two-Phase Flow
% We solve the two-phase system using a sequential splitting in which the
% pressure and fluxes are computed by solving the flow equation and then
% held fixed as the saturation is advanced according to the transport
% equation. This procedure is repeated for a given number of time steps
% (here we use 60 equally spaced time steps amounting to 50 % PV Injected).
% The error introduced by this splitting of flow and transport can be
% reduced by iterating each time step until e.g., the residual is below a
% certain user-prescribed threshold (this is not done herein).

pv     = poreVolume(G,G.rock);
nt     = 60;
Time   = 0.5*(sum(pv)/state.wellSol(1).flux);
dT     = Time/nt;
dTplot = Time/3;
N      = fix(Time/dTplot);

pvi = zeros(nt,1);
sol = cell(nt,1);

t  = 0;
count = 1; title('Saturation');
colorbar off; colormap(flipud(winter));
while t < Time,
    state = implicitTransport(state, G, dT, G.rock, fluid, 'wells', W, 'Trans', T,'verbose',true);
    
    % Check for inconsistent saturations
    assert(max(state.s(:,1)) < 1+eps && min(state.s(:,1)) > -eps);

    % Plot saturation
    delete(hp)
    hp = plotCellData(G,state.s,state.s>0); drawnow; pause(.1);
    
    % Update solution of pressure equation.
    state  = incompTPFA(state, G, T, fluid, 'wells', W, 'use_trans',true);
    sol{count,1} = state;
    
    t = t + dT;
    pvi(count) = 100*(sum(state.wellSol(1).flux)*t)/sum(pv);
    count = count + 1;
    
end

%% Plot saturations
clf, plotToolbar(G,cellfun(@(x) struct('s', x.s(:,1)), sol,'UniformOutput', false));
line(fl(:,1:2:3)',fl(:,2:2:4)',1e-3*ones(2,size(fl,1)),'Color','r','LineWidth',0.5);
hold on
plot3(xw(1,1),xw(1,2),1e-3,'xk','LineWidth',2,'MarkerSize',8);
plot3(xw(2,1),xw(2,2),1e-3,'ok','LineWidth',2,'MarkerSize',8);
hold off
colormap(flipud(winter)); caxis([0 1]); axis equal tight;

% <html>
% <p><font size="-1">
% Copyright 2009-2016 TU Delft and SINTEF ICT, Applied Mathematics.
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
