%% Basic Cartesian fracture grid
% This script illustrates the basic use of the DFM module.
% We create a simple network of two fractures in a Cartesian 2D grid and
% add the fracuters to the grid with a hybrid representation.
% Transmissibilities are then computed with a two-point flux expression,
% and solve a tracer transport problem.
%
% It is assumed that the path is set so that functions in the DFM module
% will be chosen instead of the core MRST functions. If this is not taken
% care of, strange error messages will probably result

mrstModule add dfm incomp
gravity reset off;

% Initialize a grid
Nx = [10 10];
G = cartGrid(Nx,[10 10]);
G = computeGeometry(G);

%   Create a fracture at y = 1, and between x = 1 and x = 8
fracFaces = find(all([G.faces.centroids(:,2) == 1 , G.faces.centroids(:,1) >= 1 , ...
    G.faces.centroids(:,1) < 8],2));

%   Create a fracture at x = 8, and between y = 1 and y = 10
fracFaces = [ fracFaces ; find(all([G.faces.centroids(:,1) == 8 , ...
    G.faces.centroids(:,2) >= 1 , G.faces.centroids(:,2) < 10],2))];

%   Mark the face as a fracture face
G.faces.tags=zeros(G.faces.num,1);
G.faces.tags(fracFaces) = 1;

%   Assign aperture
apt = zeros(G.faces.num,1);
aperture = 0.001;
apt(fracFaces) = aperture;

%   Add the hybrid cells
G = addhybrid(G,G.faces.tags > 0,apt);

% Plot the grid and the fracture
figure
plotGrid_DFM(G)
plotFractures(G)
axis equal, axis off

%% Add physical parameters

% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);

nCells = G.cells.num;

% Define permeability and porosity
rock = makeRock(G, milli*darcy, 0.01);

% The fracture permeability is computed from the parallel plate assumption,
% which states that the permeability is aperture squared divided by 12.
rock.perm(hybridInd,:) = aperture^2/12;
rock.poro(hybridInd) = 0.5;

% Create fluid object. The fluids have equal properties, so this ammounts
% to the injection of a tracer
fluid = initSimpleFluid('mu' , [   1,  1]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   1]);

% Compute TPFA transmissibilities
T = computeTrans_DFM(G,rock,'hybrid',true);

% Transmissibilites for fracture-fracture connections are computed in a
% separate file
[G,T2] = computeHybridTrans(G,T);

wellRate = 1;
% Injection well in lower left corner, production in upper right
W = addWell([],G,rock,1,'type','rate','val',wellRate,'comp_i',[1 0],'InnerProduct','ip_tpf');
W = addWell(W,G,rock,prod(Nx),'type','rate','val',-wellRate,'comp_i',[0 1],'InnerProduct','ip_tpf');

state = initState(G,W,0,[0 1]);


%% Solve a single-phase problem, compute time of flight, and plot
state = incompTPFA_DFM(state,G,T,fluid,'wells',W,'c2cTrans',T2);

figure
plotCellData_DFM(G,state.pressure)
plotFractures(G,hybridInd,state.pressure)
axis equal, axis off

%% Then solve tracer transport

t = 0;

% End of simulation
endTime = sum(0.1 * poreVolume(G,rock) / wellRate);

% Since the two fluids have equal properties, the pressure solution is time
% independent, and the transport equation can be solved for the entire
% simulation time at once. For visualization purposes, we split the
% interval anyhow
numSteps = 5;
dt = endTime / numSteps;

iter = 1;

while t < endTime
    state = explicitTransport_DFM(state,G,t + dt,rock,fluid,'wells',W);

    iter = iter + 1;
    t = t + dt;

    clf
    plotCellData_DFM(G,state.s(:,1));
    plotFractures(G,hybridInd,state.s(:,1));
    axis equal, axis off
    title(['Water saturation at t = ' num2str(t) 's']);
    colorbar

    pause(.1)
end

%% Copyright notice
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
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
