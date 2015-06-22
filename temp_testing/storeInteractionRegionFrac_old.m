function [CG,CGf] = storeInteractionRegionFrac(CG,CGf,CGm,varargin)
% storeInteractionRegionFrac defines and stores connectivity or grid
% topology based internal and external (matrix) interaction regions for a
% fracture coarse grid
%
% SYNOPSIS:
%   [CG,CGf] = storeInteractionRegionFrac(CG, CGf, CGm)
%   [CG,CGf] = storeInteractionRegionFrac(CG, CGf, CGm, 'pn1', pv1, ...)
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
%           1 - Classical: The matrix support will include all local
%               supports that exclude the corresponding fracture coarse
%               cells. In other words, it will be all the interaction
%               regions that each matrix cell penetrated by the fracture
%               belongs to.
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
%   sysMatrix - Fine scale system (transmissibility) matrix A (see
%               incompTPFA). Must be specified if use_metis or
%               partition_frac are true.
%
%   levels    - levels of connectivity based agggregation desired if 'type'
%               is equal to 4.
%
%   theta     - User defined tolerance for defining acceptability of a
%               coupling (connection) strength for aggregation.
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
%   storeInteractionRegion, getRsbGrids_HFM, storeFractureInternalInteractionRegion

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
    'type'            , 4                                            , ...
    'levels'          , max(ceil(0.2/CG.parent.fractureAperture),5)  , ...
    'theta'           , 0.25                                         , ...
    'excludeBoundary' , true                                         , ...
    'fullyCoupled'    , true );

opts = merge_options(opts, varargin{:});
G = CG.parent; Gm = G.Matrix; p = CG.partition;

dof_matrix = CGm.cells.num;
%
A = getConnectivityMatrix(getNeighbourship(G, 'Topological'));
CGf = storeFractureInternalInteractionRegion(CGf, G, A,...
    'simpleInteractionRegion',true,'callStoreInteractionRegion',G.griddim==3);
%
CG.cells.centers = [CGm.cells.centers; CGf.cells.centers+G.Matrix.cells.num];
CG.cells.interaction = cell(numel(CG.cells.num),1);
CG.cells.interaction(1:dof_matrix,1) = CGm.cells.interaction;
frac_cells = G.Matrix.cells.num+1:G.cells.num;
pfrac = unique(p(frac_cells));

switch opts.type

%---------------------------------Case 1----------------------------------%

    case 1
        cg_matrix = max(p(1:G.Matrix.cells.num));
        IR = zeros(G.cells.num,cg_matrix);
        for i = 1:cg_matrix
            IR(CG.cells.interaction{i,1},i) = 1;
        end
        IR = IR(1:G.Matrix.cells.num,:);
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        for i = 1:numel(pfrac)
            fcg_fcells = frac_cells(p(frac_cells)==pfrac(i)).';
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            add = []; ir_add = [];
            for j = 1:numel(matrix_cells)
                ir_add = [ir_add,find(IR(matrix_cells(j),:))]; %#ok
            end
            ir_add = unique(ir_add);
            for k = 1:numel(ir_add)
                add = [add;CG.cells.interaction{ir_add(k),1}]; %#ok
            end
            CG.cells.interaction{pfrac(i),1} = unique(add);
        end
        
%---------------------------------Case 2----------------------------------%    

    case 2
        mv = mean(CG.cells.volumes); ml = mv^(1/CG.griddim); % Change to max to cover greater area
        for i = 1:numel(pfrac)
            fcg_fcells = frac_cells(p(frac_cells)==pfrac(i)).';
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
            K = delaunayTriangulation(X);
            in = pointLocation(K, G.cells.centroids);
            interior = CGf.cells.interaction{pfrac(i)-dof_matrix,1}+G.Matrix.cells.num;
            CG.cells.interaction{pfrac(i),1} = unique([interior;find(~isnan(in))]);
        end

%---------------------------------Case 3----------------------------------%        

    case 3
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        for i = 1:numel(pfrac)
            fcg_fcells = frac_cells(p(frac_cells)==pfrac(i)).';
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            coarse_p = unique(p(matrix_cells));
            % secondary check to include any coarse cells surrounded by
            % coarse cells with a fracture in them
            nbrs = [];
            for j = 1:numel(coarse_p)
                faces = CG.cells.faces(CG.cells.facePos(coarse_p(j)):CG.cells.facePos(coarse_p(j)+1)-1,1);
                fnbr = unique(CG.faces.neighbors(faces,:));fnbr = fnbr(fnbr>0);
                fnbr = fnbr(~ismember(fnbr,coarse_p));
                nbrs = [nbrs;fnbr]; %#ok
            end
            nbrs = unique(nbrs);
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
            CG.cells.interaction{pfrac(i),1} = unique(...
                [CG.cells.interaction{pfrac(i),1}; add; ...
                find(ismember(CG.partition(1:G.Matrix.cells.num),coarse_p))]);
        end
        
%--------------------------------Case 4-----------------------------------%
    case 4
        nm = Gm.cells.num;
        for i = 1:numel(pfrac)

            pm = max(CGm.partition);
            dist = double(p == pfrac(i));
            for j = 1:opts.levels
                dist = A*dist;
            end
            CG.cells.interaction{pfrac(i),1} = find(dist>0);
        end
        for i = 1:numel(pfrac)
            start = find(p==pfrac(i)); add = start;
            for j = 1:opts.levels
                startnew = [];
                for k = 1:numel(start)
                    temp = transpose(find(A(start(k),:)));
                    temp = temp(~ismember(temp,add) & ~ismember(temp,startnew));
                    startnew = [startnew;temp]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            CG.cells.interaction{pfrac(i),1} = unique(add);
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
        temp2 = CGf.cells.interaction{j,1};
        remove2 = [CGf.cells.centers(1:j-1);CGf.cells.centers(j+1:end)];
        CGf.cells.interaction{j,1} = setdiff(temp2,remove2);
        %
        temp = CG.cells.interaction{i,1};
%         temp = [temp(temp<=Gm.cells.num); CGf.cells.interaction{j,1} + Gm.cells.num];
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
