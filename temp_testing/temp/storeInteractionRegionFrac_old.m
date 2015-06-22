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
%   including interaction in the matrix as well as within each fracture
%   network. This function internally calls storeIRFracInterior which
%   defines internal basis supports for each fracture coarse cell within
%   the corresponding fracture network.
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
%   type      - Type of fracture basis support in matrix. Must be one of
%               the following numeric values
%               1 - Classical: The matrix support will include all local
%                   supports that exclude the corresponding fracture coarse
%                   cells. In other words, it will be all the interaction
%                   regions that each matrix cell penetrated by the
%                   fracture belongs to.
% 
%               1.1 - Similar to 1 but picking only the coarse cells and
%                     not the entire interaction regions of those coarse
%                     cells. 
%
%               2 - Distance based definition of interaction region
%                   obtained by expanding the convex hull of a fracture
%                   coarse cell by the average size of a matrix coarse cell
%                   and picking all fine cells inside this new larger area.
%               
%               3 - MPFA-like interaction region picking only neighboring
%                   matrix coarse cells of the corresponding fracture
%                   coarse cell.
%
%               4 - Connectivity based definition of interaction region.
%                   Matrix connections of the internal interaction region
%                   of a fracture  coarse cell are extracted. Connections
%                   of these cells in the matrix are aggregated upto a user
%                   defined level to form the fracture basis support in the
%                   matrix.
%
%               4.1-4.4 - Combination of 4 with aggresive (Multigrid -
%                         like) local connectivity strength based
%                         aggregation. See Vanek et al, Computing, 1996 for
%                         an example.
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
% NOTE:
%   Experimental results show that all fracture interaction region types
%   computed in storeInteractionRegionFrac through the fracIRtype option
%   produce similar convergence rates and for a homogeneous multiscale
%   problem with marginal improvements in solution accuracy. However, the
%   basis functions and their iterative convergence performance may vary.
%   Graph connectivity based algorithms, more specifically fracIRtype=4
%   have proven to be simple and most efficient.
%
% SEE ALSO:
%   storeInteractionRegion, getRsbGrids_HFM, addCoarseCenterPointsFrac

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
    'sysMatrix'       , []                                           , ...
    'levels'          , max(ceil(0.2/CG.parent.fractureAperture),5) , ...
    'theta'           , 0.25                                         , ...
    'excludeBoundary' , true                                         , ...
    'fullyCoupled'    , true );

opts = merge_options(opts, varargin{:});
G = CG.parent; Gm = G.Matrix; p = CG.partition;

if opts.type >= 4
    assert(~isempty(opts.sysMatrix),'System matrix for the given grid ''G'' is not specified.');
    if opts.type>4
    if std(G.Matrix.rock.perm)==0 || (max(G.Matrix.rock.perm)/min(G.Matrix.rock.perm))<100
        warning(sprintf(['Connection strength based interaction regions ',...
            'for low-contrast systems will produce dense multiscale operators ',...
            'with bad condition numbers.\nDefaulting to aggregating upto ',...
            num2str(opts.levels),' levels of connections starting from ',...
            'matrix cells containing embedded fracture coarse cells.']));%#ok
        opts.type = 4;
    end
    end
end

dof_matrix = CGm.cells.num;
callStoreInteractionRegion = G.griddim==3;
CGf = storeIRFracInterior(CGf,G,opts.sysMatrix,'simpleInteractionRegion',true,'callStoreInteractionRegion',callStoreInteractionRegion);
CG.cells.centers = [CGm.cells.centers; CGf.cells.centers+G.Matrix.cells.num];
CG.cells.interaction = cell(numel(CG.cells.num),1);
CG.cells.interaction(1:dof_matrix,1) = CGm.cells.interaction;
frac_cells = G.Matrix.cells.num+1:G.cells.num;
cg_frac = unique(p(frac_cells));

switch opts.type
%
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
        for i = 1:numel(cg_frac)
            fcg_fcells = frac_cells(p(frac_cells)==cg_frac(i)).';
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            add = []; ir_add = [];
            for j = 1:numel(matrix_cells)
                ir_add = [ir_add,find(IR(matrix_cells(j),:))]; %#ok
            end
            ir_add = unique(ir_add);
            %-----------------------See what's happening------------------%
%             ptemp = p; ptemp(~ismember(p,ir_add))=max(ir_add)+1;
%             figure; plotGrid(G,'EdgeAlpha',0.05,'FaceColor','none'); axis tight off
%             outlineCoarseGrid(G.Matrix,ptemp(1:G.Matrix.cells.num));
%             plotGrid(G,matrix_cells,'EdgeAlpha',0,'FaceColor','b');
%             plotGrid(G,fcg_fcells,'EdgeAlpha',0,'FaceColor','y');
            %-------------------------------------------------------------%
            for k = 1:numel(ir_add)
                add = [add;CG.cells.interaction{ir_add(k),1}]; %#ok
            end
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%--------------------------------Case 1.1---------------------------------%

    case 1.1
        cg_matrix = max(p(1:G.Matrix.cells.num));
        IR = zeros(G.cells.num,cg_matrix);
        for i = 1:cg_matrix
            IR(CG.cells.interaction{i,1},i) = 1;
        end
        IR = IR(1:G.Matrix.cells.num,:);
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        for i = 1:numel(cg_frac)
            fcg_fcells = frac_cells(p(frac_cells)==cg_frac(i)).';
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            add = []; ir_add = [];
            for j = 1:numel(matrix_cells)
                ir_add = [ir_add,find(IR(matrix_cells(j),:))]; %#ok
            end
            ir_add = unique(ir_add);
            %-----------------------See what's happening------------------%
%             ptemp = p; ptemp(~ismember(p,ir_add))=max(ir_add)+1;
%             figure; plotGrid(G,'EdgeAlpha',0.05,'FaceColor','none'); axis tight off
%             outlineCoarseGrid(G.Matrix,ptemp(1:G.Matrix.cells.num));
%             plotGrid(G,matrix_cells,'EdgeAlpha',0,'FaceColor','b');
%             plotGrid(G,fcg_fcells,'EdgeAlpha',0,'FaceColor','y');
            %-------------------------------------------------------------%
            for k = 1:numel(ir_add)
                add = [add;find(p==ir_add(k))]; %#ok
            end
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%---------------------------------Case 2----------------------------------%    

    case 2
        mv = mean(CG.cells.volumes); ml = mv^(1/CG.griddim); % Change to max to cover greater area
        for i = 1:numel(cg_frac)
            fcg_fcells = frac_cells(p(frac_cells)==cg_frac(i)).';
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
            interior = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+G.Matrix.cells.num;
            CG.cells.interaction{cg_frac(i),1} = unique([interior;find(~isnan(in))]);
        %-----------------------See what's happening----------------------%
%         figure; plotGrid(G,'EdgeAlpha',0.05,'FaceColor','none'); axis tight off
%         plotGrid(G,find(~isnan(in)),'EdgeAlpha',0,'FaceColor','b');
%         plotGrid(G,fcg_fcells,'EdgeAlpha',0,'FaceColor','r');
        %-----------------------------------------------------------------%
        end
%
%---------------------------------Case 3----------------------------------%        

    case 3
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        for i = 1:numel(cg_frac)
            fcg_fcells = frac_cells(p(frac_cells)==cg_frac(i)).';
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
            %-----------------------See what's happening------------------%
%             ptemp = p; ptemp(~ismember(p,coarse_p))=max(coarse_p)+1;
%             if ~isempty(add), ptempa = p; ptempa(~ismember(p,add))=max(add)+1; end
%             figure; plotGrid(G,'EdgeAlpha',0.05,'FaceColor','none'); axis tight off
%             outlineCoarseGrid(G.Matrix,p(1:G.Matrix.cells.num),'k')
%             outlineCoarseGrid(G.Matrix,ptemp(1:G.Matrix.cells.num));
%             if ~isempty(add) outlineCoarseGrid(G.Matrix,ptempa(1:G.Matrix.cells.num),'g'); end
%             plotGrid(G,fcg_fcells,'EdgeAlpha',0,'FaceColor','b');
            %-------------------------------------------------------------%
            add = find(ismember(CG.partition(1:G.Matrix.cells.num),unique(add)));
            CG.cells.interaction{cg_frac(i),1} = unique(...
                [CG.cells.interaction{cg_frac(i),1}; add; ...
                find(ismember(CG.partition(1:G.Matrix.cells.num),coarse_p))]);
        end
%
%--------------------------------Case 4-----------------------------------%

    case 4
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        if opts.excludeBoundary
            boundary = any(G.faces.neighbors==0,2);
            facelist = 1:G.faces.num;
            bfaces = facelist( boundary);
            faces = [rldecode(1:G.cells.num,diff(G.cells.facePos),2).' G.cells.faces];
            rmcells = faces(ismember(faces(:,2),bfaces),1);
        end
        for i = 1:numel(cg_frac)
%             fcg_fcells = frac_cells(p(frac_cells)==cg_frac(i)).';
            fcg_fcells = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+nm;
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
                    temp = temp(~ismember(temp,add) & ~ismember(temp,startnew));
                    startnew = [startnew;temp]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            if opts.excludeBoundary
                add = add(~ismember(add,rmcells));
            end
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%-------------------------------Case 4.1----------------------------------%

    case 4.1
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        diagAm = zeros(nm,1);
        for i = 1:nm
            diagAm(i) = abs(sum(Am(i,:)));
        end
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        for i = 1:numel(cg_frac)
            fcg_fcells = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+nm;
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            start = [];
            for j = 1:numel(matrix_cells)
                start = [start;transpose(find(Am(matrix_cells(j),:)))]; %#ok
            end
            add = unique([fcg_fcells;matrix_cells;start]);
            start = unique(start(~ismember(start,matrix_cells)));
            for j = 1:opts.levels-1
                startnew = [];
                for k = 1:numel(start)
                    temp = transpose(find(Am(start(k),:)));
                    temp = temp(~ismember(temp,add) & ~ismember(temp,startnew));
                    startnew = [startnew;temp]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            while ~isempty(startnew)
                startnew = [];
                for j = 1:numel(start)
                    conn = transpose(find(Am(start(j),:)));
                    strong = [];
                    for k = 1:numel(conn)
                        if Am(start(j),conn(k))>opts.theta*sqrt(diagAm(start(j))*diagAm(conn(k)))
                            strong = [strong;conn(k)]; %#ok
                        end
                    end
                    strong = strong(~ismember(strong,add) & ~ismember(strong,startnew));
                    startnew = [startnew;strong]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%-------------------------------Case 4.2----------------------------------%       

    case 4.2
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        diagAm = zeros(nm,1);
        for i = 1:nm
            diagAm(i) = abs(sum(Am(i,:)));
        end
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        for i = 1:numel(cg_frac)
            fcg_fcells = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+nm;
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            start = [];
            for j = 1:numel(matrix_cells)
                start = [start;transpose(find(Am(matrix_cells(j),:)))]; %#ok
            end
            add = unique([fcg_fcells;matrix_cells;start]);
            start = unique(start(~ismember(start,matrix_cells)));
            startnew = start; 
            while ~isempty(startnew)
                startnew = [];
                for j = 1:numel(start)
                    conn = transpose(find(Am(start(j),:)));
                    strong = [];
                    for k = 1:numel(conn)
                        if Am(start(j),conn(k))>=opts.theta*sqrt(diagAm(start(j))*diagAm(conn(k)))
                            strong = [strong;conn(k)]; %#ok
                        end
                    end
                    strong = strong(~ismember(strong,add) & ~ismember(strong,startnew));
                    startnew = [startnew;strong]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
%             pir = ones(G.cells.num,1); pir(unique(add)) = 2;
%             pir = processPartition(G,pir);
%             add = find(pir==2);
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%-------------------------------Case 4.3----------------------------------%

    case 4.3
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        stencil = zeros(nm,1);
        diagAm = stencil;
        for i = 1:nm
            stencil(i) = numel(find(Am(i,:))); % stencil size without diagonal
            diagAm(i) = abs(sum(Am(i,:)));
        end
        stencil = max(stencil);
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        for i = 1:numel(cg_frac)
            fcg_fcells = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+nm;
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            start = [];
            for j = 1:numel(matrix_cells)
                start = [start;transpose(find(Am(matrix_cells(j),:)))]; %#ok
            end
            add = unique([fcg_fcells;matrix_cells;start]);
            start = unique(start(~ismember(start,matrix_cells)));
            for j = 1:opts.levels-1
                startnew = [];
                for k = 1:numel(start)
                    temp = transpose(find(Am(start(k),:)));
                    temp = temp(~ismember(temp,add) & ~ismember(temp,startnew));
                    startnew = [startnew;temp]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            strength = opts.theta*max(diagAm(unique(add)))/stencil;
            while ~isempty(startnew)
                startnew = [];
                for j = 1:numel(start)
                    strong = transpose(find(Am(start(j),:)>=strength));
                    strong = strong(~ismember(strong,add) & ~ismember(strong,startnew));
                    startnew = [startnew;strong]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end
%
%-------------------------------Case 4.4----------------------------------%       

    case 4.4
        nm = G.Matrix.cells.num;
        Am = opts.sysMatrix(1:nm,1:nm);
        Am = Am-diag(diag(Am)); % Remove diagonal entries, it also includes fracture contribution
        stencil = zeros(nm,1);
        for i = 1:nm
            stencil(i) = numel(find(Am(i,:))); % stencil size without diagonal
        end
        stencil = max(stencil);
        nnc = G.nnc.cells;
        nnc = nnc(ismember(nnc(:,1),1:G.Matrix.cells.num),:);
        Am = abs(Am);
        for i = 1:numel(cg_frac)
            fcg_fcells = CGf.cells.interaction{cg_frac(i)-dof_matrix,1}+nm;
            matrix_cells = unique(nnc(ismember(nnc(:,2),fcg_fcells),1));
            start = []; diagAm = zeros(size(matrix_cells));
            for j = 1:numel(matrix_cells)
                diagAm(j) = abs(sum(Am(matrix_cells(j),:))); 
                start = [start;transpose(find(Am(matrix_cells(j),:)))]; %#ok
            end
            strength = opts.theta*max(diagAm)/stencil; % Local definition of strength in the vicinity of the network
            add = unique([fcg_fcells;matrix_cells;start]);
            start = unique(start(~ismember(start,matrix_cells)));
            startnew = start; 
            while ~isempty(startnew)
                startnew = [];
                for j = 1:numel(start)
                    strong = transpose(find(Am(start(j),:)>strength));
                    strong = strong(~ismember(strong,add) & ~ismember(strong,startnew));
                    startnew = [startnew;strong]; %#ok
                end
                start = startnew;
                add = [add;start]; %#ok
            end
%             pir = ones(G.cells.num,1); pir(unique(add)) = 2;
%             pir = processPartition(G,pir);
%             add = find(pir==2);
            CG.cells.interaction{cg_frac(i),1} = unique(add);
        end

%-------------------------------------------------------------------------%        

    otherwise
        error('Wrong option for defining fracture interaction regions. Argument ''type'' in function call must be an integer between 1-4.');
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
        temp = [temp(temp<=Gm.cells.num); CGf.cells.interaction{j,1} + Gm.cells.num];
        remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
        CG.cells.interaction{i,1} = setdiff(setdiff(temp,remove1),rmcells);
    else
        if opts.fullyCoupled
            temp = CG.cells.interaction{i,1};
            remove1 = [CG.cells.centers(1:i-1);CG.cells.centers(i+1:end)];
            CG.cells.interaction{i,1} = setdiff(temp,remove1);
        end
    end
end
return


function CGf = storeIRFracInterior(CGf,G,A,varargin)

opt = struct(...
    'simpleInteractionRegion'    , true  , ...
    'callStoreInteractionRegion' , false );
    
opt = merge_options(opt, varargin{:});

if opt.callStoreInteractionRegion == true
    CGf = storeInteractionRegion(CGf,...
    'adjustCenters'              , false  , ...
    'skipSingleCellBlocks'       , false  , ... 
    'simpleCellGrouping'         , true   , ...
    'ensureConnected'            , true   , ...
    'largeBasis'                 , true   , ...  % true is better
    'localTriangulation'         , false  , ...  % local triangulation will not work for fractures
    'useMultipoint'              , false  );  
elseif opt.simpleInteractionRegion == true
    CGf = addCoarseCenterPointsFrac(CGf,'option','useCoarseFaceCentroids','meantype','arithmetic'); % IR centers
    Af = A(G.Matrix.cells.num+1:end,G.Matrix.cells.num+1:end);
    Af = Af - diag(diag(Af));
    p = CGf.partition;
    for i = 1:max(p)
        rmcenters = [CGf.cells.centers(1:i-1);CGf.cells.centers(i+1:end)];
        fcg_fcells = find(p==i);
        irf = fcg_fcells;
        start = irf;
        while ~isempty(start)
            start2 = [];
            for j = 1:numel(start)
                add = find(Af(start(j),:));
                add = add(~ismember(add,irf));
                add = add(~ismember(add,rmcenters));
                start2 = [start2;add']; %#ok
            end
            start = start2;
            irf = [irf;start]; %#ok
        end
        CGf.cells.interaction{i,1} = unique(irf);
    end
end
return