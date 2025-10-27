function [bc, bd_cells, bd_faces, wellCells,entire_boundary] = introduceBoundaryConditionsFromSchdule(G, states, schedule, topZ, bottomZ, caseType, data_ratio)
% Define boundary conditions based on specified case type for MRST simulations.
%
% Inputs:
%   G         - MRST grid structure.
%   states    - Array of states for each simulation timestep.
%   schedule  - MRST schedule structure containing control settings.
%   topZ      - Z-coordinate of the top faces for boundary conditions.
%   bottomZ   - Z-coordinate of the bottom faces for boundary conditions.
%   caseType  - String specifying boundary condition case:
%               'onlywells'          - Apply control only on well faces.
%               'wellsPlusSeismic'    - Apply control on well faces and random faces from top/bottom.
%               'wellPlusSeismicTop'  - Apply control on well faces and random faces from top only.
%   data_ratio - Number of random faces selected from top/bottom faces for additional control.
%
% Outputs:
%   bc        - Boundary condition structure for MRST.
%   bd_cells  - Selected boundary cells based on specified conditions.
%   bd_faces  - Selected boundary faces based on specified conditions.

% Initialization
bc = cell(length(schedule.step.val), 1);
W = schedule.control.W;

% Identify top and bottom faces based on Z-coordinates
if G.griddim > 2 &&all(G.cells.centroids(:, 3)~= G.cells.centroids(:, 3))
    topFaces = find(abs(G.faces.centroids(:, 3) - topZ) < eps);
    bottomFaces = find(abs(G.faces.centroids(:, 3) - bottomZ) < eps);
    entire_boundary = boundaryFaces(G);
else
    topFaces = find(abs(G.faces.centroids(:, 2) - topZ) < eps);
    bottomFaces = find(abs(G.faces.centroids(:, 2) - bottomZ) < eps);
        leftFaces = find(abs(G.faces.centroids(:, 1) - topZ) < eps);
   rightFaces = find(abs(G.faces.centroids(:, 1) - bottomZ) < eps);
   entire_boundary = [topFaces; bottomFaces];
end
allBoundaryFaces = [topFaces; bottomFaces];

% Gather all cells in wells
wellCells = unique(vertcat(W.cells));

% Identify well boundary faces and cells
[wellbdFaces, wellbdCells] = findWellBoundaryFaces(G,wellCells, allBoundaryFaces);
entire_boundary = setdiff(entire_boundary, wellbdFaces);
seismicFaces =[];
seismicCells=[];
% Select boundary control case
switch caseType
    case 'onlywells'
        bd_cells = wellbdCells;
        bd_faces = wellbdFaces;

        % Apply boundary conditions
        for i = 1:length(schedule.step.val)
            bc{i} = addBC([], bd_faces, 'pressure', states{i}.pressure(bd_cells), 'sat', states{i}.s(bd_cells, :));
        end
        
    case 'wellsPlusSeismic'
        % Randomly select a subset of top/bottom faces
        numDataFaces = numel(allBoundaryFaces);
        assert(data_ratio <= numDataFaces, 'Data ratio exceeds available boundary faces.');
        randIndices = randperm(numDataFaces, data_ratio);
        seismicFaces = allBoundaryFaces(randIndices);
        seismicFaces = setdiff(seismicFaces, wellbdFaces);
        seismicCells = unique(sum(G.faces.neighbors(seismicFaces, :), 2));
        bd_faces = unique([wellbdFaces; seismicFaces]);
        bd_cells = sum(G.faces.neighbors(bd_faces, :), 2);
        % Apply boundary conditions
        for i = 1:length(schedule.step.val)
            bc{i} = addBC([], seismicFaces, 'pressure', states{i}.pressure(seismicCells), 'sat', states{i}.s(seismicCells, :));
        end
        bd_cells = seismicCells;
        bd_faces = seismicFaces;
        
    case 'wellPlusSeismicTop'
        % Randomly select a subset of top faces only
        numTopFaces = numel(topFaces);
        assert(data_ratio <= numTopFaces, 'Data ratio exceeds available top faces.');
        randIndices = randperm(numTopFaces, data_ratio);
        seismicFaces = topFaces(randIndices);
        seismicFaces = setdiff(seismicFaces, wellbdFaces);
        seismicCells = unique(sum(G.faces.neighbors(seismicFaces, :), 2));
        bd_faces = [wellbdFaces; seismicFaces];
        bd_cells = sum(G.faces.neighbors(bd_faces, :), 2);
        % Apply boundary conditions
        for i = 1:length(schedule.step.val)
            bc{i} = addBC([], bd_faces, 'pressure', states{i}.pressure(bd_faces), 'sat', states{i}.s(bd_faces, :));
        end
    otherwise
        warning('No valid boundary control case specified; no boundary control applied.');
        bd_cells = [];
        bd_faces = [];
        return;
end
end

function [wellBoundaryFaces, wellBoundaryCells] = findWellBoundaryFaces(G, wellCells, allBoundaryFaces)
% FINDTOPBOTTOMWELLBOUNDARYFACES Find boundary faces of well cells at the top and bottom
% of the MRST grid based on the provided bottomTopFaces vector.
%
% Inputs:
%   G              - MRST grid structure.
%   wellCells      - Vector of cell indices considered as well cells.
%   bottomTopFaces - Vector of face indices that correspond to the top and bottom of the grid.
%
% Outputs:
%   wellBoundaryFaces - Indices of boundary faces belonging to well cells at top and bottom.
%   wellBoundaryCells - Indices of well cells associated with these boundary faces.

% Get neighbors for the specified bottomTopFaces
neighbors = sum(G.faces.neighbors(allBoundaryFaces, :),2);

% Create a logical mask for faces connected to well cells
isWellFace = ismember(neighbors, wellCells);

% Select boundary faces connected to well cells
wellBoundaryFaces = allBoundaryFaces(isWellFace);
% treat bottom faces separately
maxDepth = max(G.cells.centroids(:, 3));
wellCentroids = G.cells.centroids(wellCells, :);

if all(max(wellCentroids(:,3))<maxDepth)
    % Find associated well cells

    % Step 1: Extract x and y coordinates of well cells
    wellX = wellCentroids(:, 1);  % x-coordinates
    wellY = wellCentroids(:, 2);  % y-coordinates

    % Step 2: Find bottom/top face centroids
    faceCentroids = G.faces.centroids(allBoundaryFaces, :);

    % Step 3: Identify faces below wells (same x, y and lower z)
    underWellBoundaryMask = false(size(allBoundaryFaces));  % Initialize mask
if G.griddim > 2

    for i = 1:numel(wellX)
        % Find faces with matching (x, y) and z-coordinate below the current well
        sameXY = (abs(faceCentroids(:, 1) - wellX(i)) < eps) & ...
            (abs(faceCentroids(:, 2) - wellY(i)) < eps);
        belowWell = faceCentroids(:, 3) > wellCentroids(i, 3);

        % Mark faces that are directly below the well in the mask
        underWellBoundaryMask = underWellBoundaryMask | (sameXY & belowWell);
    end
else
     for i = 1:numel(wellX)
        % Find faces with matching (x, y) and z-coordinate below the current well
        sameX =  abs(faceCentroids(:, 2) - wellY(i)) < eps;
        belowWell = faceCentroids(:, 2) < wellCentroids(i, 2);

        % Mark faces that are directly below the well in the mask
        underWellBoundaryMask = underWellBoundaryMask | (sameX & belowWell);
     end
end

    % Step 4: Select the faces directly below the wells
    underWellBoundaryFaces = allBoundaryFaces(underWellBoundaryMask);
    wellBoundaryFaces = unique([wellBoundaryFaces;underWellBoundaryFaces]);
end
wellBoundaryCells = sum(G.faces.neighbors(wellBoundaryFaces, :),2);
end

