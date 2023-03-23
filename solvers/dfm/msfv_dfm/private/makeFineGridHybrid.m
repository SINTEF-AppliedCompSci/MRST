function g_fine = makeFineGridHybrid(g_dual,f_vertices,f_fractures,dx,dy,box,varargin)
% make a fine grid where the coarse dual grid and the fine-scale fractures
% are constrains in the grid. Hybrid cells can be added using
%
% SYNOPSIS:
% g = addHybridCells(g,h);
%
% PARAMETERS:
%   g_dual: a mrst-grid structure acting as constrains for the finescale
%           model
%   f_fractures: fine-scale fractures (with tags)
%   f_vertices:  the vertices defining the fine-scale fractures
%   dx,dy:       distance between points
%   box:         bounding box of the domain
% OPTIONAL:
%   precision:   the precision of the fine-scale grid
%
% Copyright 2011-2012 University of Bergen, 2013 IRIS AS
%
% This file is licensed under the GNU General Public License v3.0.
opt = struct ('precision', 1e-8 ...
    );
warning('off','makeFineGridHybrid:Option:Unsupported')
opt = merge_options(opt, varargin{:});

% merge the coarse constraints with the fine-scale fractures
vertices = g_dual.nodes.coords;
edges = [double(reshape(g_dual.faces.nodes,2,[])') g_dual.faces.tags (1:g_dual.faces.num)'];
if ~isempty(f_fractures)
    f_fractures = [f_fractures zeros(size(f_fractures,1),1)];
    [vertices,edges] = mergeFractures(vertices,edges,f_vertices,f_fractures,box,opt.precision);
end

% remove duplicated edges.
[vertices,edges] = removeDuplicates(vertices,edges);

% create the points for the triangulation
[X,Y] = meshgrid(box(1,1):dx:box(2,1),box(1,2):dy:box(2,2));
Y(1:end-1,2:2:end) = Y(1:end-1,2:2:end)+0.5 * dy;
p     = [X(:), Y(:)];

% take away extra points that are too close to the line to be of any use
tol = sqrt(3)/6 * dx;
p = remove_closepoints (vertices, edges, p, tol);


% we need to take care of the crossings our self in order to
% keep the tags
[vertices, edges] = removeFractureIntersections (vertices, edges, box, opt);

% distribute points along the edges.
[vertices,edges,~] = partition_edges (vertices,edges,dx,box,opt);

% add up the vertices
vertices = [vertices ; p];

% triangulate and return a mrst grid
g_fine = triangulate(vertices,edges(:,1:2),edges(:,3:end));

% tags are stored as a matrix the meaning of the columns are
% specified in tagmark.
% for instance g_fine.faces.tags(:,g_fine.faces.tagmark.frac)
% returns the fracture tags.
g_fine.faces.tagmark.frac = 1;
g_fine.faces.tagmark.boundary = 2;
g_fine.faces.tagmark.map2c = 3;
