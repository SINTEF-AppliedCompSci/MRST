function [CG,CGf] = storeFractureInteractionRegion(CG,CGf,CGm,varargin)
% storeFractureInteractionRegion defines and stores connectivity or grid
% topology based internal and external (matrix) support regions for a
% given matrix and fracture coarse grid.
%
% SYNOPSIS:
%   [CG,CGf] = storeFractureInteractionRegion(CG, CGf, CGm)
%   [CG,CGf] = storeFractureInteractionRegion(CG, CGf, CGm, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function computes the support for fracture basis functions
%   including interaction in the matrix as well as internal to each
%   fracture network.
%
% REQUIRED PARAMETERS:
%
%   CG  - Coarsened global grid 'G' containing fracture and matrix cells.
%
%   CGf - Coarsened fracture grid 'Gf' defined by assembleFracGrid.
%
%   CGm - Matrix coarse grid.
%
% OPTIONAL PARAMETERS:
%
%   levels    - levels of connectivity used to define the support region
%               for a fracture coarse block. In other words, the
%               topological radius of the fracture support region.
%
%   fullyCoupled - Coupled fracture and matrix basis functions. Allows
%                  matrix support to extend into fractures. Increases
%                  accuracy but reduces speed.
%
%   excludeBoundary - Excludes fine cells on the boundary from support
%                     regions of internal coarse blocks. Improves accuracy,
%                     especially when boundary conditions are specified.
%
%
%%%%%%%%%%%%%%%%%%%-----------EXPERIMENTAL-----------%%%%%%%%%%%%%%%%%%%%%%
%
%
%   removeCenters    - Logical value to essentially enforce a basis
%                      function value of 1 at each coarse node. May improve
%                      convergence rate in an iterative multiscale
%                      strategy.
%
%   coarseNodeOption - Possible methods for optimizing coarse node location
%                      inside fractures.
%
% RETURNS:
%   CG - Coarse grid for combined fracture and matrix grid 'G' with basis
%        supports and coarse node indices in added fields
%        'G.cells.interaction' and 'G.cells.centers respectively'. 
%
%   CGf - Fracture coarse grid with similar added information as CG but
%         only for the internal fracture grid Gf as returned by
%         assembleFracGrid.
%
% SEE ALSO:
%   storeInteractionRegion, getRsbGridsHFM, 
%   storeFractureInternalInteractionRegion, addCoarseCenterPointsFrac

%{
Copyright 2009-2015: TU Delft and SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opts = struct('levels'          , 8      , ...
              'excludeBoundary' , true   , ...
              'fullyCoupled'    , true   , ...
              'removeCenters'   , true   , ...
              'coarseNodeOption', 'useFineCellCentroids')  ;

opts = merge_options(opts, varargin{:});
G = CG.parent; Gm = G.Matrix; p = CG.partition;
dof_matrix = CGm.cells.num;
%
A = getConnectivityMatrix(getNeighbourship(G, 'Topological'));
%
CGf = addCoarseCenterPointsFrac(CGf,'option',opts.coarseNodeOption);
CG.cells.centers = [CGm.cells.centers; CGf.cells.centers+G.Matrix.cells.num];
CG.cells.interaction = cell(CG.cells.num,1);
CG.cells.interaction(1:dof_matrix,1) = CGm.cells.interaction;
pfrac = unique(p(Gm.cells.num+1:G.cells.num));

%-------------------------------------------------------------------------%        

hwb = waitbar(0,'Define and store fracture connectivity..');
for i = 1:numel(pfrac)
    F = gridCellFaces(CG,pfrac(i));
    N = CG.faces.neighbors(F,:);
    N = N(N~=0);
    N = N(N~=pfrac(i));
    centers = CG.cells.centers(N);
    if G.griddim>2
        dist = double(p == pfrac(i));
    else
        CGf = storeFractureInternalInteractionRegion(CGf, G, A,...
            'simpleInteractionRegion', true,...
            'callStoreInteractionRegion', G.griddim==3);
        dist = zeros(size(p));
        dist(CGf.cells.interaction{pfrac(i)-dof_matrix,1} + Gm.cells.num) = 1;
    end
    flag = 0;
    for j = 1:opts.levels
        dist = A*dist;
        if G.griddim>2
            if all(ismember(centers,find(dist>0)))
                flag = 1;
                ir = find(dist>0);
                inFrac = ir(ir>Gm.cells.num);
                dist = A*dist;
                ir = find(dist>0);
                inMat = ir(ir<=Gm.cells.num);
                ir = [inFrac;inMat];
                break;
            end
        end
    end
    if flag == 0, ir = find(dist>0); end
    CG.cells.interaction{pfrac(i),1} = ir;
    waitbar(i/numel(pfrac),hwb);
end
close(hwb);

%-------------------------------------------------------------------------%        

if opts.excludeBoundary
    boundary = any(Gm.faces.neighbors==0,2);
    facelist = 1:Gm.faces.num;
    bfaces = facelist(boundary);
    faces = [rldecode(1:Gm.cells.num,diff(Gm.cells.facePos),2).' Gm.cells.faces];
    rmcells = faces(ismember(faces(:,2),bfaces),1);
end

hwb = waitbar(0,'Define and store fracture support..');
for i = 1:CG.cells.num
    if i>dof_matrix
        j = i-dof_matrix;
        % CGf only
        if isfield(CGf.cells, 'interaction') && opts.removeCenters
            temp2 = CGf.cells.interaction{j,1};
            remove2 = [CGf.cells.centers(1:j-1);CGf.cells.centers(j+1:end)];
            CGf.cells.interaction{j,1} = setdiff(temp2,remove2);
        end
        % CG (fracture part)
        temp = CG.cells.interaction{i,1};
        if G.griddim == 2
            if ~isfield(CGf.cells, 'interaction')
                CGf = storeFractureInternalInteractionRegion(CGf, G, A,...
                      'simpleInteractionRegion',true,...
                  	  'callStoreInteractionRegion', G.griddim==3);
            end
            temp = [ temp(temp<=Gm.cells.num); ...
                     CGf.cells.interaction{j,1} + Gm.cells.num];
        end
        if opts.removeCenters
            remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
            CG.cells.interaction{i,1} = setdiff(temp,remove1);
        else
            CG.cells.interaction{i,1} = temp;
        end
        if opts.excludeBoundary
            CG.cells.interaction{i,1} = setdiff(CG.cells.interaction{i,1}, rmcells);
        end
    else
        % CG (matrix part)
        if opts.fullyCoupled && opts.removeCenters
            temp = CG.cells.interaction{i,1};
            remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
            CG.cells.interaction{i,1} = setdiff(temp,remove1);
        end
    end
    waitbar(i/CG.cells.num,hwb);
end
close(hwb);
return
