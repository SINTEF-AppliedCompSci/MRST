function data = setupReservoirBoundaryConditions(model, states, schedule, bcType, seismicFaces)
%SETUPRESERVOIRBOUNDARYCONDITIONS Configure boundary conditions for reservoir simulation
%
% SYNOPSIS:
%   data = setupReservoirBoundaryConditions(model, states, schedule, bcType, seismicFaces)
%
% PARAMETERS:
%   model            - Reservoir model
%   states           - Cell array of state vectors for each timestep
%   schedule         - Simulation schedule
%   bcType           - Boundary condition type:
%                     'wells'    - Well boundaries only
%                     'top'      - Top surface boundaries
%                     'seismic'  - User-provided seismic faces
%                     'combined' - Wells + top + seismic
%   seismicFaces     - (Optional) List of seismic face indices
%
% RETURNS:
%   data             - Structure containing all boundary condition data:
%                     .bc            - Boundary condition structure
%                     .bdCells      - Boundary cell indices
%                     .bdFaces      - Boundary face indices
%                     .wellCells     - Well cell indices
%                     .AllFaces - All boundary faces in grid

% Input validation
G = model.G;
validateattributes(G, {'struct'}, {}, mfilename, 'G', 1);
validateattributes(states, {'cell'}, {}, mfilename, 'states', 2);
validateattributes(schedule, {'struct'}, {}, mfilename, 'schedule', 3);
validTypes = {'wells', 'top', 'seismic', 'combined'};
if ~any(strcmpi(bcType, validTypes))
    error('Invalid boundary condition type. Choose from: wells, top, seismic, combined');
end

if strcmpi(bcType, 'seismic') && nargin < 5
    error('Seismic faces must be provided for ''seismic'' boundary type');
end

% Initialize outputs
data.bc = cell(numel(states), 1);
data.AllFaces = findBoundaryFaces(G);
data.wellCells = getWellCells(schedule);

% Find well-connected boundary faces
[wellFaces, wellCellsOnBoundary] = findWellBoundaryFaces(G, data.wellCells);
wellCells = data.wellCells;
% Setup boundary conditions based on type
switch lower(bcType)
    case 'wells'
        % Well boundaries only
        data.bdFaces = wellFaces;
        data.bdCells = wellCellsOnBoundary;

    case 'top'
        % Top surface boundaries
        topFaces = identifyBcFaces(G, 'Upper', [], [], []);
        data.bdFaces = topFaces;
        data.bdCells = sum(G.faces.neighbors(topFaces, :), 2);

    case 'seismic'
        % User-provided seismic faces
        validateattributes(seismicFaces, {'numeric'}, {'vector', 'positive', '<=', G.faces.num}, ...
            mfilename, 'seismicFaces', 5);
        neighbors = G.faces.neighbors(seismicFaces, :);
        isWellFace = ismember(neighbors(:,1), wellCells) | ismember(neighbors(:,2), wellCells);
        seismicFaces = seismicFaces(~isWellFace);
        data.bdFaces = seismicFaces;
        data.bdCells = sum(G.faces.neighbors(seismicFaces, :), 2);

    case 'combined'
        % Combination of wells, top, and seismic faces
        topFaces = identifyBcFaces(G, 'Upper', [], [], []);
        if nargin < 5
            seismicFaces = [];
        end
        data.bdFaces = unique([wellFaces; topFaces; seismicFaces(:)]);
        data.bdCells = sum(G.faces.neighbors(data.bdFaces, :), 2);
end

% Create boundary conditions for each timestep
for i = 1:numel(states)
    if any(model.gravity)
        gRhoDz = computeGravity(model, states{i});
        facePressure = states{i}.pressure(data.bdCells) + gRhoDz{2}(data.bdCells);
    else
        facePressure = states{i}.pressure(data.bdCells);
    end
    data.bc{i} = addBC([], data.bdFaces, 'pressure', facePressure, ...
        'sat', states{i}.s(data.bdCells, :));
end

end

%% Helper functions (remain exactly the same as in your original code)
function faces = findBoundaryFaces(G)
% Find all boundary faces in grid
faces = find(any(G.faces.neighbors == 0, 2));
end

function cells = getWellCells(schedule)
% Get all well cells from schedule
cells = [];
for c = 1:numel(schedule.control)
    if isfield(schedule.control(c), 'W')
        cells = [cells; vertcat(schedule.control(c).W.cells)];
    end
end
cells = unique(cells);
end

function [wellFaces, wellCells] = findWellBoundaryFaces(G, wellCells)
% Find boundary faces connected to well cells
topFaces = identifyBcFaces(G, 'Upper', [], [], []);
neighbors = G.faces.neighbors(topFaces, :);
isWellFace = ismember(neighbors(:,1), wellCells) | ismember(neighbors(:,2), wellCells);
wellFaces = topFaces(isWellFace);
wellCells = sum(G.faces.neighbors(wellFaces, :), 2);
end

function bcfaces = identifyBcFaces(G, side, i1, i2, i3)
bcfaces = boundaryFaceIndices(G,side,i1,i2,i3);
end


function gRhoDz = computeGravity(model, state)
% to compute gravity-adjusted face pressures for each phase and state.
%
% INPUT:
%   model  - Reservoir model (e.g., ThreePhaseBlackOilModel).
%   states - Cell array of state objects (one per time step).
%   prop   - Struct with optional settings:
%              * weight: Scaling factor for gravity terms (default: 1).
%              * saturationWeighting: If true, weight density by phase saturation (default: false).
%
% OUTPUT:
%   gRhoDz - Cell array of size (numel(states), nph) containing gravity terms for each phase.

if nargin < 3
    prop = struct('weight', 1, 'saturationWeighting', false);
end

act = model.getActivePhases();
nph = sum(act);
gRhoDz = cell(1, nph);
nf = size(model.operators.N, 1);

if norm(model.gravity) > 0 && nf > 0
    gdz = -model.getGravityGradient();
    nm = model.getPhaseNames();
    rho = state.PVTProps.Density;
    rho = expandMatrixToCell(rho);
    avg = model.operators.faceAvg;

    for i = 1:nph
        if prop.saturationWeighting
            s = model.getProp(state, ['s', nm(i)]);
            rhof = avg(s .* rho{i}) ./ max(avg(s), 1e-8);
        else
            rhof = avg(rho{i});
        end
        gRhoDz{i} = rhof .* gdz;
    end
else
    [gRhoDz{:}] = deal(zeros(nf, 1));
end

% Apply weighting if provided
if isfield(prop, 'weight') && ~isempty(prop.weight)
    for i = 1:nph
        gRhoDz{i} = prop.weight * gRhoDz{i};
    end
end

end