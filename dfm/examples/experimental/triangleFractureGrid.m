%% Unstructured fracture grid created using the Triangle software
% Triangulate a fractured domain using the Triangle software. The fractures
% are read from a odp-file (LibreOffice presentation).
%
% Read the fracture network from the file simpleSystem.odp.
% Some fractures have the same color, these will be assigned the same tag
% if the color option is passed to build_fractures_mod. Note that the
% input color values (specified in RGB values) must be the same as those
% used in the odp-file; the filter is not robust with respect to inexact
% color specifications. Sorry..
% If colors are not specified, each fracture gets assigned a unique tag.

mrstModule add dfm incomp


fractures = fullfile(mrstPath('dfm'), 'simpleSystem.odp');

if exist(fractures, 'file') ~= 2,
   error(['Fracture file ''simpleSystem.odp'' ', ...
          'does not exist on this system']);
end

fracInfo = struct('fractures', fractures, ...
    'color',[0 1 0;1 0 0], ...
    'precision', 1e-2);

% build the fractures
[vertices,fractures,box] = build_fractures_mod(fracInfo);

%% Add boundaries and create grid using Triangle

% Add the boundary to the point cloud passed to Triangle.
pts = [vertices ; ...
    box(1,1) box(1,2) ;...
    box(2,1) box(1,2); ...
    box(2,1) box(2,2);...
    box(1,1) box(2,2)];

numVert = size(vertices,1);

boundaryEdges = [ 1 2 ; 2 3 ; 3 4 ; 4 1] + numVert;

% Assign negative edge markers to the boundaries
boundaryEdges = [boundaryEdges , (-1 : -1 : -4)'];

edges = [fractures ; boundaryEdges];

% Triangulate using Triangle
mrstModule add triangle
G = mex_triangleGrid(pts, edges(:,1:2), 'maxArea', 0.3,'segmentmarkers',edges(:,3));
G = computeGeometry(G);

%% Recover tag information
% For some reason (quite possibly user stupidity, but issues with the mex
% compiler is also a candidate), Triangle does not provide reasonable edge
% tags. Recover the tags by identifying faces that have both nodes on the
% fracture.
% If the tags from Triangle works, ie the content in G.faces.tag is
% sensible, you can simply replace the next lines with
%
%   G.faces.tags = G.faces.tag
%
% and move on.
%
% And yes, tags are used in a sub-optimal way in this module..

% Recover faceNodes
faceNodes = reshape(G.faces.nodes,2,[])';

G.faces.tags = zeros(G.faces.num,1);

% Loop over all fractures, compare faceNode indices with points that are on
% the fracture
% Note that the reason this works is that Triangle adds points instead of
% moving the fractures to get a reasonable grid.
for iter = 1 : size(fractures,1)

    % Endpoints of all fractures. We can find points on the boundary in a
    % similar matter.
    x0 = vertices(fractures(iter,1),1); y0 = vertices(fractures(iter,1),2);
    x1 = vertices(fractures(iter,2),1); y1 = vertices(fractures(iter,2),2);

    % Compute distance from nodes to this fracture, find those that are on the
    % fracture (within rounding errors)
    dist = distance_to_closest_line(G.nodes.coords(:,1),G.nodes.coords(:,2),x0,y0,x1,y1);
    onFrac = find(dist < sqrt(eps));

    % Find faces with both point on the fracture
    faceOnThisLine = all(ismember(faceNodes,onFrac),2);
    G.faces.tags(faceOnThisLine) = fractures(iter,3);
end

%% Create hybrid cells in fratures

% Apertures of all faces. Zero aperture for non-fractures
aperture = zeros(G.faces.num,1);

% Fractures of family 1 are larger both in length and aperture than the
aperture(G.faces.tags == 1) = 10^-3;
aperture(G.faces.tags == 2) = 10^-4;


% Add hybrid cells along the fractures
G = addhybrid(G,G.faces.tags > 0,aperture);

% Plot the grid and the fracture. When Triangle is used for gridding, the
% fractures are represented as straight lines.
figure
plotGrid_DFM(G)
plotFractures(G)
axis equal, axis off

%% Set parameters

% Find indices of hybrid cells
hybridInd = find(G.cells.hybrid);
nCells = G.cells.num;

% Define permeability and porosity
rock.perm = milli * darcy * ones(nCells,2);
rock.poro = 0.01 * ones(nCells,1);

% Much higher values are used in the fracture
rock.perm(hybridInd,:) = aperture(G.cells.tags(hybridInd)).^2/12 * [1 1];
rock.poro(hybridInd) = 1;

% Create fluid object. Non-linear rel perms, resident fluid has the higher
% viscosity
fluid = initSimpleFluid('mu' , [   2,  2]*centi*poise     , ...
    'rho', [1014, 859]*kilogram/meter^3, ...
    'n'  , [   1,   10]);

%% Define wells

% Well rate
wellRate = 1;

% Place both injection and production wells close to the fracture network
[~,injInd] = min(sum([G.cells.centroids(:,1) - 7 , G.cells.centroids(:,2) - 7.5].^2,2));
[~,prodInd] = min(sum([G.cells.centroids(:,1) - 9 , G.cells.centroids(:,2) - .2].^2,2));

W = addWell([],G,rock,injInd,'type','rate','val',wellRate,'comp_i',[1 0]);
W = addWell(W,G,rock,prodInd,'type','rate','val',-wellRate,'comp_i',[0 1]);

state = initState(G,W,0,[0 1]);

%% Compute transmissibilities with both MPFA and TPFA

% Transmissibilities with TPFA
T_TPFA = computeTrans_DFM(G,rock,'hybrid',true);
[G,T_TPFA_2] = computeHybridTrans(G,T_TPFA);
state_TPFA = incompTPFA_DFM(state,G,T_TPFA,fluid,'wells',W,'c2cTrans',T_TPFA_2);

% Transmissibilities with MPFA
T_MPFA = computeMultiPointTrans_DFM(G,rock,'hybrid',true);
[G,T_MPFA] = computeHybridMPTrans(G,T_MPFA);
state_MPFA = incompMPFA_DFM(state,G,T_MPFA,fluid,'wells',W,'cellConnections',true);



%% Solve two-phase problem

% Comments are sparse here, for explanations see examples in the core files
% of MRST.
% The plotting is designed to highlight differences between two- and
% multi-point fluxes

t = 0;

% Simulation endtime
T = 0.03 * sum(poreVolume(G,rock)) / wellRate;

numTimeSteps = 30;
dt = T / numTimeSteps;

numSatPlots = 3;
dtSatPlot = T / numSatPlots;

nextSatPlot = dtSatPlot;

% Store concentration in the outlet
wellConcentration = zeros(numTimeSteps + 1,2);


while t < T
    % Transport solve with both methods. We use implicit transport to avoid
    % severe time step restrictions due to high flow rates and small cells
    state_TPFA =  implicitTransport_DFM(state_TPFA,G,t + dt,rock,fluid,'wells',W);
    state_MPFA =  implicitTransport_DFM(state_MPFA,G,t + dt,rock,fluid,'wells',W);

    wellConcentration(iter + 1,1) = state_TPFA.s(prodInd,1);
    wellConcentration(iter + 1,2) = state_MPFA.s(prodInd,1);

    % Update pressure solution
    state_TPFA = incompTPFA_DFM(state_TPFA,G,T_TPFA,fluid,'wells',W,'c2cTrans',T_TPFA_2);
    state_MPFA = incompMPFA_DFM(state_MPFA,G,T_MPFA,fluid,'wells',W,'cellConnections',true);

    iter = iter + 1;
    t = t + dt;

    % Plotting
    if t >= nextSatPlot
        figure,
        subplot(2,1,1)
        plotCellData_DFM(G,state_TPFA.s(:,1));
        plotFractures(G,hybridInd,state_TPFA.s(:,1));
        title('TPFA')
        axis off, axis equal


        subplot(2,1,2)
        plotCellData_DFM(G,state_MPFA.s(:,1));
        plotFractures(G,hybridInd,state_MPFA.s(:,1));
        title('MPFA')
        axis off, axis equal



        nextSatPlot = nextSatPlot + dtSatPlot;
    end
end

% Finally plot water cuts in the producer
figure, hold on
h = plot(wellConcentration(:,1));
h = [h plot(wellConcentration(:,2),'r')];
legend(h,'TPFA','MPFA','FontSize',16)

