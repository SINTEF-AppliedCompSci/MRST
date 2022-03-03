function model = CreateGrid_3D(model)
    % 3-D grid in x-direction
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
    number_of_y_cells = 10; %the real number will be this - 1
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
        y = linspace(0, sqrt(area), number_of_y_cells); model.grid.y = y;
        z = linspace(0, sqrt(area), number_of_y_cells); model.grid.z = z;
        dxs = [dxLeft,dx * ones(1,nCells),dxRight]; model.grid.dxs = dxs;
        dys = dy * ones(nCells+2,1); model.grid.dys = dys;
        dzs = dz * ones(nCells+2,1); model.grid.dxz = dzs;        
    else
        x = (0:dx:nCells*dx)'; model.grid.x = x;
        y = linspace(0, sqrt(area), number_of_y_cells); model.grid.y = y;
        z = linspace(0, sqrt(area), number_of_y_cells); model.grid.z = z;
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
    
    
    % find index of inlet face
    inlet_mask = G.cells.centroids(:,1) == G.cells.centroids(1,1);
    outlet_mask = G.cells.centroids(:,1) == G.cells.centroids(G.cartDims(1,1), 1);
    inlet_frontcells_mask = G.cells.centroids(:,1) == G.cells.centroids(2,1);
    outlet_backcells_mask = G.cells.centroids(:,1) == G.cells.centroids(G.cartDims(1,1) - 1, 1);

    grid.G.inlet_mask = inlet_mask;
    grid.G.outlet_mask = outlet_mask;
    grid.G.inlet_frontcells_mask = inlet_frontcells_mask;
    grid.G.outlet_backcells_mask = outlet_backcells_mask;
    
    % create satnum regions
    if (isfield(simulation,'bCells'))
        satNum       = ones(grid.G.cells.num, 1);
        satNum(inlet_mask)    = 2; satNum(outlet_mask) = 2;
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