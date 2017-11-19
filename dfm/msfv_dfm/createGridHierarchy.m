function [g_dual,g_fine]= createGridHierarchy(coarseFractures,fineFractures,varargin)
%CREATEGRIDHIERARCHY Generates the fine-scale grid and the dual-coarse grid
%
% SYNOPSIS
% [g_dual,g_fine,box,opt] = createGridHierarchy(coarseGrids,fineFractures,varargin)
%
% PARAMETERS
%
%   - coarseFractures - a structure defining the dual coarse grid
%                   The structure is an input into build_fractures_mod
%                   which reads fractures drawn a openoffice file (.odp)
%                   or given in a .mat file
%
%   - fineFractures - a cell of structures defining the fine fractures
%                   The structures are input into build_fractures_mod (see above)
%   OPTIONAL
%   - load          - path to load grid (if empty nothing is loaded)
%   - save          - path to save grid (if empty nothing is saved)
%   - numElements   - number of fine-scale cells (the final number of cells
%                     are usally 10% more then this.
%   - apertures     - a vector with apertures. One for each tag. (if empty the
%                     the apertures specified in the input structures are
%                     used.
%   - shrinkfactor  - The width of the edge cells that are not fractures are
%                     computed such that these edge cells have approximatly
%                     the same size as the inner cells. In order to use a
%                     hybrid representation of these celle their width must
%                     be much smaller than the neighboring cells. A
%                     shrinkfactor can be used to adjust the width of these
%                     cells.
%
%   NB Varargin is pasted directly into subrutines, and options can thus be
%   pasted directly to its subrutines.
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.

% default options
opt = struct ('load', '', ...
    'save', '', ...
    'numElements', 100000, ...
    'apertures', [], ...
    'shrinkfactor', 10, ...
    'box', [], ...
    'precision', 1e-3, ...
    'merge', true ...
    );
warning('off','createGridHiearchy:Option:Unsupported')
opt = merge_options(opt, varargin{:});

% load a pre made grid
if ~isempty(opt.load)
    load(opt.load)
    return
end

% in no aperture is specified we find them in the input structures
if isempty(opt.apertures)

    for i = 1: numel(coarseFractures)
        apertures(coarseFractures{i}.tag) = coarseFractures{i}.aperture;
    end

    for i = 1: numel(fineFractures)
        apertures(fineFractures{i}.tag) = fineFractures{i}.aperture;
    end
    opt.apertures = apertures(:);
    clear apertures
end

% make the coarse dual Grid
disp('Make coarse grid')

% get the large-scale fractures.
[vertices,fractures] = readFractures(coarseFractures,opt.box,opt.precision);

% make the coarse dual grid
g_dual = makeCoarseDualGrid(vertices,fractures,opt.box,opt.merge,opt.precision);

% make the fine grid
disp('Make fine grid')

% ideal spacing in x and y directions
[opt.dx,opt.dy] = getIdealSpacing(opt.box,opt.numElements);

% There should be enough cells to resolve the coarse fractures
% If error either adjust the precision of the coarse fractures
% or add more fine-scale cells.
assert(min(g_dual.cells.volumes)>opt.dx*opt.dy)

% get the fine-scale fractures.
[vertices,fractures] = readFractures(fineFractures,opt.box,opt.precision);

% make the fine grid
[g_fine] = makeFineGridHybrid(g_dual,vertices,fractures,opt.dx,opt.dy,opt.box,varargin{:});

% create hybrid cells for the fine fractures.
g_fine = addHybridCells(g_fine,g_dual,opt);

% save the grid.
if ~isempty(opt.save)
    save(opt.save,'g_dual','g_fine','coarseFractures','fineFractures','opt')
end


end

function g_fine = addHybridCells(g_fine,g_dual,opt)
% Add hybrid cells for the edge cells, node cells and
% fine-scale fractures.

%% Add hybrid cells for the edge cells and node cells

% map to coarse grid
ctags = g_fine.faces.tags(:,g_fine.faces.tagmark.map2c);

% boundary tags
btags = g_fine.faces.tags(:,g_fine.faces.tagmark.boundary);

% fracture tags
ftags = g_fine.faces.tags(:,g_fine.faces.tagmark.frac);

% make hybrid cells for all coarse constraints
% that are not at the boundary
% and all fine-scale fractures
hybrid = ctags>0 & btags == 0 | (ftags>0 & ctags == 0);

%% find the width
hf = zeros(g_fine.faces.num,1);

% compute the width of the non-fracture edge cells

% unit normals
normals = bsxfun(@times,g_dual.faces.normals,1./g_dual.faces.areas);

% the width of the face cells (we divide by two since it is rectangles not
% triangles. In this way the areas should be similar. A shrink factor is
% used to make sure the width are smaller then the neighboring cells
hc = sqrt(sum((normals.*(ones(g_dual.faces.num,1)*[opt.dx,opt.dy])).^2,2))./opt.shrinkfactor;

% associate the coarse width to the fine
hf(ctags>0 & btags == 0) = hc(ctags(ctags>0 & btags == 0));

% the width for the fracture edges are just their aperture

% we map the apertures to the tags
[utags,~,map] = unique(g_fine.faces.tags(:,g_fine.faces.tagmark.frac));
numfrac = sum(utags>0);
numshift = numel(utags)-numfrac;
utags = utags(utags>0);
a_tmp = [zeros(numshift,1);opt.apertures(utags)];
apertures = a_tmp(map);
hf(ftags>0) = apertures(ftags>0);

% and stores them in the structure for later usage
g_fine.faces.apertures = apertures;

%% create hybrid cells of kind 1 and kind 2(node cells)
g_fine = addhybrid(g_fine,hybrid,hf,'addHybrid2Cells',true);

end

function [dx,dy] = getIdealSpacing(box,numElements)
% get the spacing of nodes given a a the number of elements

% find the area per. element.
x_size = box(2)-box(1);
y_size = box(4)-box(3);
total_area = x_size * y_size;
area_per_elem = total_area / numElements;


% initially we setup the grid so that every two rows are moved to the
% middle between the points on the row above. this makes the distance
% between the projection of the point above to this row and the point
% itself 1/2*x. the hypotenuse is x, since we want equilateral triangles
% which makes the last side in the right-angle triangle (which fills
% half the element) sqrt(3)/2*x. the area of an element thus relates to
% x such that it is 1/2 * x * sqrt(3)/2*x
dx = sqrt (4 * area_per_elem / sqrt (3));
dy = sqrt (3) / 2 * dx;

% adjust so that we get a round nice number of points in each direction
nx = max (1, round (x_size / dx));
dx = x_size / nx;
ny = max (1, round (y_size / dy));
dy = y_size / ny;
end

