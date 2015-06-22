function [CG,CGf] = storeFractureInteractionRegion(CG,CGf,CGm,varargin)
% storeInteractionRegionFrac defines and stores connectivity or grid
% topology based internal and external (matrix) interaction regions for a
% fracture coarse grid
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
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   type - Type of fracture basis support in matrix. Must be one of
%          the following numeric values
%           1 - Classical: Fracture support in matrix will include all
%                          coarse blocks belonging to all matrix
%                          interaction regions a fracture coarse cell is
%                          part of.
% 
%           2 - Distance based definition of interaction region
%               obtained by expanding the convex hull of a fracture coarse
%               cell by the average size of a matrix coarse cell and
%               picking all fine cells inside this new larger area.
% 
%           3 - MPFA-like interaction region picking only neighboring
%               matrix coarse cells of the corresponding fracture coarse
%               cell.
% 
%           4 - Connectivity based definition of interaction region.
%               Matrix connections upto a certain topological distance are
%               selected starting from fracture-matrix NNC connections.
%
%   levels    - levels of connectivity based agggregation desired if 'type'
%               is equal to 4.
%
%   excludeBoundary - exclude matrix cells at the boundary from fracture
%                     interaction regions. Useful when BC is applied to
%                     matrix
%
%   fullyCoupled - same as getRsbGrids_HFM
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
%   storeInteractionRegion, getRsbGrids_HFM, 
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


opts = struct(...
    'type'            , 4                                             , ...
    'levels'          , max(ceil(0.2/CG.parent.fractureAperture),5)  , ...
    'excludeBoundary' , true                                          , ...
    'fullyCoupled'    , true );

opts = merge_options(opts, varargin{:});
G = CG.parent; Gm = G.Matrix; p = CG.partition;
dof_matrix = CGm.cells.num;
%
A = getConnectivityMatrix(getNeighbourship(G, 'Topological'));
%
CGf = addCoarseCenterPointsFrac(CGf,'option','useCoarseCellCentroids','meantype','arithmetic');
CG.cells.centers = [CGm.cells.centers; CGf.cells.centers+G.Matrix.cells.num];
CG.cells.interaction = cell(CG.cells.num,1);
CG.cells.interaction(1:dof_matrix,1) = CGm.cells.interaction;
pfrac = unique(p(Gm.cells.num+1:G.cells.num));

switch opts.type

%---------------------------------Case 1----------------------------------%

    case 1
        CGf = storeFractureInternalInteractionRegion(CGf, G, A,...
            'simpleInteractionRegion',true,...
            'callStoreInteractionRegion', G.griddim==3);
        cg_matrix = max(p(1:Gm.cells.num));
        IR = zeros(G.cells.num,cg_matrix);
        for i = 1:cg_matrix
            IR(CG.cells.interaction{i,1},i) = 1;
        end
        IR = IR(1:Gm.cells.num,:);
        for i = 1:numel(pfrac)
            dist = double(p == pfrac(i)); dist = A*dist;
            mcells = find(dist>0); mcells = mcells(mcells<=Gm.cells.num);
            touch = zeros(1,Gm.cells.num); touch(mcells) = 1;
            ir_add = find(touch*IR>0);
            add = CGf.cells.interaction{pfrac(i) - dof_matrix,1} + Gm.cells.num;
            CG.cells.interaction{pfrac(i),1} = [add;find(ismember(p,ir_add))];
        end
        
%---------------------------------Case 2----------------------------------%    

    case 2
        mv = mean(CG.cells.volumes); ml = mv^(1/CG.griddim); % Change to max to cover greater area
        for i = 1:numel(pfrac)
            fcg_fcells = find(p==pfrac(i));
            faces = gridCellFaces(G,fcg_fcells);
            xcent = G.faces.centroids(faces,1); ycent = G.faces.centroids(faces,2);
            if CG.griddim == 3
                zcent = G.faces.centroids(faces,3);
                xcent = [xcent;xcent-ml;xcent+ml;xcent;xcent;xcent;xcent]; %#ok
                ycent = [ycent;ycent;ycent;ycent-ml;ycent+ml;ycent;ycent]; %#ok
                zcent = [zcent;zcent;zcent;zcent;zcent;zcent+ml;zcent-ml]; %#ok
                X = [xcent,ycent,zcent];
            else
                xcent = [xcent;xcent-ml;xcent+ml;xcent;xcent]; %#ok
                ycent = [ycent;ycent;ycent;ycent-ml;ycent+ml]; %#ok
                X = [xcent,ycent];
            end
            K = delaunayTriangulation(unique(X,'rows'));
            in = pointLocation(K, G.cells.centroids);
            CG.cells.interaction{pfrac(i),1} = find(~isnan(in));
        end

%---------------------------------Case 3----------------------------------%        

    case 3
        CGf = storeFractureInternalInteractionRegion(CGf, G, A,...
            'simpleInteractionRegion',true,...
            'callStoreInteractionRegion', G.griddim==3);
        for i = 1:numel(pfrac)
            dist = double(p == pfrac(i)); dist = A*dist;
            mcells = find(dist>0); mcells = mcells(mcells<=Gm.cells.num);
            coarse_p = unique(p(mcells));
            % secondary check to include any coarse cells surrounded by
            % coarse cells with a fracture in them
            nbrs = [];
            for j = 1:numel(coarse_p)
                faces = CG.cells.faces(CG.cells.facePos(coarse_p(j)):CG.cells.facePos(coarse_p(j)+1)-1,1);
                fnbr = unique(CG.faces.neighbors(faces,:));fnbr = fnbr(fnbr>0);
                fnbr = fnbr(~ismember(fnbr,[coarse_p;nbrs]));
                nbrs = [nbrs;fnbr]; %#ok
            end
            add = [];
            for j = 1:numel(nbrs);
                faces = CG.cells.faces(CG.cells.facePos(nbrs(j)):CG.cells.facePos(nbrs(j)+1)-1,1);
                fnbr = unique(CG.faces.neighbors(faces,:)); fnbr = fnbr(fnbr>0);
                isit = ismember(fnbr,coarse_p);
                if numel(find(isit))==numel(fnbr) || numel(find(isit))==(numel(fnbr)-1)
                    add = [add;nbrs(j)]; %#ok
                end
            end
            add = add(add>0);
            add = find(ismember(CG.partition(1:G.Matrix.cells.num),unique(add)));
            add_interior = CGf.cells.interaction{pfrac(i) - dof_matrix,1} + Gm.cells.num;
            CG.cells.interaction{pfrac(i),1} = unique(...
                [CG.cells.interaction{pfrac(i),1}; add; add_interior; ...
                find(ismember(CG.partition(1:G.Matrix.cells.num),coarse_p))]);
        end
        
%--------------------------------Case 4-----------------------------------%

    case 4
        for i = 1:numel(pfrac)
            F = gridCellFaces(CG,pfrac(i));
            N = CG.faces.neighbors(F,:);
            N = N(N~=0);
            N = N(N~=pfrac(i));
            centers = CG.cells.centers(N);
            dist = double(p == pfrac(i));
            flag = 0;
            for j = 1:opts.levels
                dist = A*dist;
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
            if flag == 0, ir = find(d>0); end
            CG.cells.interaction{pfrac(i),1} = ir;
        end
        
%-------------------------------------------------------------------------%        

    otherwise
        error('Invalid argument. Value of key ''type'' must be an integer between 1-4.');
end

if opts.excludeBoundary
    boundary = any(Gm.faces.neighbors==0,2);
    facelist = 1:Gm.faces.num;
    bfaces = facelist(boundary);
    faces = [rldecode(1:Gm.cells.num,diff(Gm.cells.facePos),2).' Gm.cells.faces];
    rmcells = faces(ismember(faces(:,2),bfaces),1);
end

for i = 1:CG.cells.num
    if i>dof_matrix
        j = i-dof_matrix;
        %
        if isfield(CGf.cells, 'interaction')
            temp2 = CGf.cells.interaction{j,1};
            remove2 = [CGf.cells.centers(1:j-1);CGf.cells.centers(j+1:end)];
            CGf.cells.interaction{j,1} = setdiff(temp2,remove2);
        end
        %
        temp = CG.cells.interaction{i,1};
        if G.griddim == 2
            CGf = storeFractureInternalInteractionRegion(CGf, G, A);
            temp = [ temp(temp<=Gm.cells.num); ...
                     CGf.cells.interaction{j,1} + Gm.cells.num];
        end
        remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
        if opts.excludeBoundary
            CG.cells.interaction{i,1} = setdiff(setdiff(temp,remove1),rmcells);
        else
            CG.cells.interaction{i,1} = setdiff(temp,remove1);
        end
    else
        if opts.fullyCoupled
            temp = CG.cells.interaction{i,1};
            remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
            CG.cells.interaction{i,1} = setdiff(temp,remove1);
        end
    end
end
return

%{
    case 4.1
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        for i = 1:numel(pfrac)
            fcg_fcells = CGf.cells.interaction{pfrac(i)-dof_matrix,1}+nm;
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            start = []; diagAm = zeros(size(matrix_cells));
            for j = 1:numel(matrix_cells)
                diagAm(j) = abs(sum(Am(matrix_cells(j),:))); 
                start = [start;transpose(find(Am(matrix_cells(j),:)))]; %#ok
            end
            add = unique([fcg_fcells;matrix_cells;start]);
            start = unique(start(~ismember(start,matrix_cells)));
            for j = 1:opts.levels-1
                startnew = [];
                for k = 1:numel(start)
                    temp = transpose(find(Am(start(k),:)));
                    if any(ismember(temp,CG.cells.centers)),continue, end
                    temp = temp(~ismember(temp,add) & ~ismember(temp,startnew));
                    startnew = [startnew;temp]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            CG.cells.interaction{pfrac(i),1} = unique(add);
        end
%}
