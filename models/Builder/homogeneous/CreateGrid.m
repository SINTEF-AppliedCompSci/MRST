function model = CreateGrid(model)
% <keywords>
%
% Purpose : create the simulation grid
%
% Syntax :
%   model = CreateGrid(model)
%
% Input Parameters :
%   model: struct output from the configure function
%
% Return Parameters :
%   model: struct cointaing the grid properties
%
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
    % 1-D grid in x-direction
    experiment = model.experiment;
    simulation = model.simulation;
    length = experiment.geometry.length.value;
    diameter = experiment.geometry.diameter.value;
    nCells = simulation.nCells.value;
    dx = length / nCells;
    area = pi * diameter ^ 2 / 4;
    model.grid.area = area;
    dy = sqrt(area); dz = dy;
    model.grid.dx = dx; model.grid.dy = dy; model.grid.dz = dz;
    if (isfield(simulation,'bCells'))        
        bCells = simulation.bCells.value;
        x = zeros(nCells+3,1);
        dxLeft  = bCells * dx;
        dxRight = bCells * dx;
        x(1) = 0;
        x(2) = dxLeft; 
        x(3:end-1) = dxLeft + (dx:dx:nCells*dx)';
        x(end) = x(end-1) + dxRight;
        model.grid.x = x;
        y = [0;dy]; model.grid.y = y;
        z = [0;dz]; model.grid.z = z;
        dxs = [dxLeft,dx * ones(1,nCells),dxRight]; model.grid.dxs = dxs;
        dys = dy * ones(nCells+2,1); model.grid.dys = dys;
        dzs = dz * ones(nCells+2,1); model.grid.dxz = dzs;        
    else
        x = (0:dx:nCells*dx)'; model.grid.x = x;
        y = [0;dy]; model.grid.y = y;
        z = [0;dz]; model.grid.z = z;
        dxs = dx * ones(nCells,1); model.grid.dxs = dxs;
        dys = dy * ones(nCells,1); model.grid.dys = dys;
        dzs = dz * ones(nCells,1); model.grid.dzs = dzs;
    end
   
    % create MRST grid
%     G = cartGrid([nCells+2,1,1], [max(x),max(y),max(z)]*meter);    
    G = tensorGrid(x, y, z);
    G = computeGeometry(G);
    G = gridAddHelpers(G);
    grid.G = G;  
    
    % create satnum regions
    if (isfield(simulation,'bCells'))
        satNum       = ones(grid.G.cells.num, 1) * 2;
        satNum(1)    = 1; satNum(end) = 1;
        grid.satNum  = satNum;     
    end
    
    % display the grid
    %---------------------------------------
%     figTitle   = 'Grid';
%     xLabel     = 'x';
%     yLabel     = 'y'; 
%     zLabel     = 'z';
%     figStyle   = 'docked';
%     figTag     = figTitle;
%     fig = figure('Name',        figTitle, ...
%                  'Tag',         figTag,   ...
%                  'NumberTitle', 'off',     ...
%                  'WindowStyle', figStyle );
%     figure(fig);
%     xlabel(xLabel); ylabel(yLabel); zlabel(zLabel);
%     xLim = [min(x) max(x)]; yLim = [min(y) max(y)]; zLim = [min(z) max(z)];
%     xlim(xLim); ylim(yLim); zlim(zLim);
%     plotGrid(G,'FaceColor',[.7 .7 1]); view(3); axis equal;
%     title('Grid Representation');
    %---------------------------------------
    
    model.grid = grid;
end