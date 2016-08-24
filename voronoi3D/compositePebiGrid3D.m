function [G,F] = compositePebiGrid3D(celldim, pdims, varargin)
% Construct a 2D composite Pebi grid. A cartesian background grid that is 
% refined around faults and wells.
%
% SYNOPSIS:
%   G = compositePebiGrid(resGridSize, pdims)
%   G = compositePebiGrid(...,'Name1',Value1,'Name2',Value2,...)
%
% PARAMETERS:
%   celldim           - [xSize,ySize, zSize] Size of the reservoir grid 
%                       cells in x, y, and z direction.
%   pdim              - Vector, length 2, [xmax, ymax, zmax], of physical
%                       size in units of meters of the computational domain 
%
%   faultSurf         - OPTINAL
%                       Default value empty. A nfx1 cell array of 
%                       triangulations. Each triangulation should be a 
%                       surface embedded in 3D. The grid will conform to 
%                       the surfaces.
%   faultRho         - OPTINAL.
%                       Default value sqrt(sum(pdims)). faultRho is a cell
%                       array of length nf. Each element in faultRho is a 
%                       function. The function must map from a nx3 array to
%                       nx1 array. faultRho(faultSurf{i}.Points) define the
%                       radius of the balls placed along fault i.
%
% RETURNS:
%   G                - Valid grid definition.  
%
% EXAMPLE:
%   
%
% SEE ALSO:
%   compositePebiGrid, pebi, createFaultGridPoints, createWellGridPoints.

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}  

% Set options
opt = struct('faultSurf', {}, ...
             'faultRho', {});         

opt = merge_options(opt, varargin{:});

% Create fault sites
F = createFaultGridPoints3D(opt.faultSurf,opt.faultRho);

% Create reservoir sites
x = pdim(1); y = pdim(2); z = pdim(3);

xa = linspace(0, x, celldim(1)+1);
ya = linspace(0, y, celldim(2)+1);
za = linspace(0, z, celldim(3)+1);

bdr   = [ 0, 0, 0;  ...
          x, 0, 0;  ...
          x, y, 0;  ...
          0, y, 0;  ...
          0, 0, z;  ...
          x, 0, z;  ...
          x, y, z;  ...
          0, y, z];

        [X,Y,Z] = ndgrid(xa,ya,za);
rSites = [X(:), Y(:), Z(:)];

% Remove conflict points
nR = size(rSites,1);
rGs = zeros(nR,1);
rPri = zeros(nR,1);

[rSites, removed] = faultSufCond3D(rSites,F.c.CC,F.c.R);
rGs = rGs(~removed);
rPri = rPri(~removed);

% Create grid
pts = [F.f.pts;rSites];
gs  = [F.f.Gs;rGs];
pri = [F.f.pri;rPri];
%[pts,removed] = removeConflictPoints3(pts,4*gs, pri );
%pts = [X(:),Y(:),Z(:)];
bdrDT = delaunayTriangulation(bdr);
Gd = compositePebiGrid3D(pts,bdrDT);

end


