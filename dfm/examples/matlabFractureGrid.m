%% Unstructured fracture grid created by Matlab functions
% A somewhat advanced example, where a fracture network is gridded
% with triangles using a fairly simple approach. The flow equations are
% then discretized with both two- and multi-point flux expressions, and a
% two-phase transport problem is solved.

mrstModule add dfm incomp
% Create an initial distribution of points
N     = 20;
N1    = N;
N2    = N;
[X, Y] = meshgrid(0:1:N1, 0:1:N2);
X     = sqrt(3) / 2 * X;
Y(:,2:2:end) = Y(:,2:2:end)+0.5;
p     = [X(:), Y(:)];

% Estimate of the grid size
h = .3;

% Define endpoints of the fractures.
% If a nertex here is also contained in p, one of them should be removed to
% avoid mapping errors in the subsequent Delaunay triangulation
vertices = [ 3 1 ; 16 10 ;...
    3 15.7 ; 16 18 ;...
    15 4 ; 7 19.5  ;...
    4 5  ; 8 13  ;...
    14 15; 13 10];

% Give the constraints wrt numbering of the vertices
constraints = reshape(1 : size(vertices,1),2,[])' ;

% Assign tags to the fractures
constraints = [constraints, (1 : size(constraints,1))'];

% We create the grid by combining the intial points with points distributed
% along the fractures, and gridding the total point cloud with a Delaunay
% algorithm constrained by the fractures. Some tricks are needed to ensure
% the resulting grid has reasonable quality, these are briedly explained
% below. Also note that this strategy does not create high-quality grids,
% for that more advanced methods are needed; see for instance another
% example in the DFM module for a grid created by triangle.

% Bounding box of the domain. Needed for legacy reasons (when the point
% that eventually defines the grid is preconditioned to improve the grid
% quality, points are added to the cloud, and the box is needed to ensure
% the extra points are within the domain).
box = [0 0; N1 N2];

% Cells close the fractures should be significantly larger than the
% apertures for the hybrid DFM approximation to be reasonable. Moreover,
% points that are too close to the fractures will also produce weired
% cells. For these reasons, remove points that are close to the fractures.
% NOTE: The distance from the fractures where we remove points should be
% significantly larger than the tolerance for the fracture location defined
% next.
perimiter = h ;
p = remove_closepoints(vertices,constraints,p,perimiter);

% The quality of the grid is improved if the fractures are alowed to
% deviate somewhat from straight lines. This is of course an approximation,
% but it is in agreement with a concept of uncertainty in the fracture
% description. This parameter is more important for more advanced
% gridding algorithms than the present one. Still, reducing this value will
% make the fractures closer to straight lines
args = struct('precision',0.01);

% Add vertices in the intersection between fractures
[vertices, constraints] = removeFractureIntersections(vertices,constraints,box,args);

numOrdPt = size(p,1);

% Concatenate the prescribed points and the fracture endpoints
p  = [p ; vertices];

% Let constraints refer to the new point distribution
constraints(:,1:2) = constraints(:,1:2) + numOrdPt;

% Add points along the fractures to get reasonable grids there as well
[p, constraints, map] = partition_edges(p,constraints,0.5,box,args);

% Bookkeeping
tags = constraints(:,3);
constraints =  constraints(:,1:2);

% Now the points are defined. Create a constrained Delauny triangulation
% on the point cloud, and create the grid
delTri     = DelaunayTri(p, constraints);
G     = triangleGrid(delTri.X, delTri.Triangulation);
G = computeGeometry(G);

%% Introduce hybrid representation of fractures

G.faces.tags = zeros(G.faces.num,1);

% Recover faceNodes, and sort the nodes columnwise (for the use of ismember
% later)
faceNodes = sort(reshape(G.faces.nodes,2,[])',2);

% Use constraints from the Delaunay triangulation, in case these have been
% changed (it shouldn't happen, and will give a warning)
% Rowwise sort to prepare for using ismember
constraints = sort(delTri.Constraints,2);

% Loop over all constraints (the segments of fractures) and find the faces
% that spans the constraint.
for iter = 1 : size(constraints)
    % ismember works here since both inputs have been sorted
    fracFace = find(ismember(faceNodes,constraints(iter,:),'rows'));

    % Assign a tag to the newly found fracture face
    G.faces.tags(fracFace) = tags(iter);
end

% Apertures of all faces. Zero aperture for non-fractures
aperture = zeros(G.faces.num,1);

% Set different apertures for each
for iter = unique(tags)'
    hit = G.faces.tags == iter;
    aperture(hit) = 10^-3 * rand(1); % This choice is random in more than one way..
end

% Add hybrid cells along the fractures
G = addhybrid(G,G.faces.tags > 0,aperture);

% Plot the grid and the fracture. Note that the fracuters are not straight
% lines.
figure
plotGrid_DFM(G)
plotFractures(G)
axis equal, axis off

%% Set parameters

% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);
nCells = G.cells.num;

% Define permeability and porosity
rock = makeRock(G, [1, 1]*milli*darcy, 0.01);

% Much higher values are used in the fracture
rock.perm(hybridInd,:) = aperture(G.cells.tags(hybridInd)).^2/12 * [1 1];
rock.poro(hybridInd) = 0.5;

% Create fluid object. Non-linear rel perms, resident fluid has the higher
% viscosity
fluid = initSimpleFluid('mu' , [   2,  2]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   10]);

%% Define wells

% Well rate
wellRate = 1;

% Place both injection and production wells close to the fracture network
injInd = delTri.pointLocation(3, 1);
prodInd = delTri.pointLocation(16, 17);

W = addWell([],G,rock,injInd,'type','rate','val',wellRate,'comp_i',[1 0],'InnerProduct','ip_simple');
W = addWell(W,G,rock,prodInd,'type','rate','val',-wellRate,'comp_i',[0 1],'InnerProduct','ip_simple');

state = initState(G,W,0,[0 1]);

%% Compute transmissibilities with both MPFA and TPFA

% Note that there is a difference in how connections between hybrid cells
% are handled between TPFA and MPFA: Both methods discretize the
% connections in files that are separate from compute(MP)Trans, but TPFA
% gives a separate list of transmissibilities for cell-to-cell connections,
% whereas MPFA modifies the standard transmissibility matrix.

% Transmissibilities with TPFA
T_TPFA = computeTrans_DFM(G,rock,'hybrid',true);
[G,T_TPFA_2] = computeHybridTrans(G,T_TPFA);
state_TPFA = incompTPFA_DFM(state,G,T_TPFA,fluid,'wells',W,'c2cTrans',T_TPFA_2);

figure
plotCellData_DFM(G,state_TPFA.pressure)
plotFractures(G,hybridInd,state_TPFA.pressure)
title('TPFA')
axis off, axis equal
% 
% % Transmissibilities with MPFA
% T_MPFA = computeMultiPointTrans_DFM(G,rock,'hybrid',true);
% [G,T_MPFA] = computeHybridMPTrans(G,T_MPFA);
% state_MPFA = incompMPFA_DFM(state,G,T_MPFA,fluid,'wells',W,'cellConnections',true);
% 
% figure
% plotCellData_DFM(G,state_MPFA.pressure)
% plotFractures(G,hybridInd,state_MPFA.pressure)
% title('TPFA')
% axis off, axis equal

%% Solve two-phase problem.

% Comments are sparse here, for explanations see examples in the core files
% of MRST.
% The plotting is designed to highlight differences between two- and
% multi-point fluxes

t = 0;

% Simulation endtime
T = 0.02 * sum(poreVolume(G,rock)) / wellRate;

numTimeSteps = 100;
dt = T / numTimeSteps;

numSatPlots = 3;
dtSatPlot = T / numSatPlots;

nextSatPlot = dtSatPlot;

% Store concentration in the outlet
wellConcentration = zeros(numTimeSteps + 1,2);

% Transport loop.
while t < T
    % Transport solve with both methods. We use implicit transport to avoid
    % severe time step restrictions due to high flow rates and small cells
    state_TPFA =  implicitTransport_DFM(state_TPFA,G,t + dt,rock,fluid,'wells',W);
%     state_MPFA =  implicitTransport_DFM(state_MPFA,G,t + dt,rock,fluid,'wells',W);

    wellConcentration(iter + 1,1) = state_TPFA.s(prodInd,1);
%     wellConcentration(iter + 1,2) = state_MPFA.s(prodInd,1);

    % Update pressure solution
    state_TPFA = incompTPFA_DFM(state_TPFA,G,T_TPFA,fluid,'wells',W,'c2cTrans',T_TPFA_2);
%     state_MPFA = incompMPFA_DFM(state_MPFA,G,T_MPFA,fluid,'wells',W,'cellConnections',true);

    iter = iter + 1;
    t = t + dt;

    % Plotting
%     if t >= nextSatPlot
        figure(1); clf
%         subplot(2,1,1)
        plotCellData_DFM(G,state_TPFA.s(:,1));
        plotFractures(G,hybridInd,state_TPFA.s(:,1));
        title('TPFA')
        axis off, axis equal
        drawnow
% 
%         subplot(2,1,2)
%         plotCellData_DFM(G,state_MPFA.s(:,1));
%         plotFractures(G,hybridInd,state_MPFA.s(:,1));
%         title('MPFA')
%         axis off, axis equal
% 
%         nextSatPlot = nextSatPlot + dtSatPlot;
%     end

end

% Finally plot water cuts in the producer
figure, hold on
plot(wellConcentration(:,1))
plot(wellConcentration(:,2),'r')
title('Producer water cut')

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


