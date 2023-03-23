function [F,fracture] = assembleFracNodes2D(G,fracture,varargin)
% assembleFracNodes2D partitions fracture lines to generate nodes for
% gridding
%
% SYNOPSIS:
%   [F,fracture] = assembleFracNodes2D(G, fracture);
%   [F,fracture] = assembleFracNodes2D(G, fracture, 'pn1', pv1, ...);
%
% DESCRIPTION:
%   This function supports 3 possible options to divide a fracture for
%   gridding. The resulting output can be thought of as a 1D grid for each
%   fracture line. From preliminary experiments, option/type 1 has proved
%   to be the cheapest way to ensure relatively uniform gridding with
%   similar resolution throughout all fractures. For each assemblyType
%   fracture intersections will be a node on each of the intersecting line.
%   This is because transmissibility at these intersections is determined
%   using the star-delta transformation (see SPE-88812-PA, Karimi-Fard et
%   al, 2004). Hence, the intersection must be a fine-scale face when the
%   intersecting fractures are represented as narrow 2D cells in the global
%   grid.
%
% REQUIRED PARAMETERS:
%
%   G        - Matrix grid structure.
%
%   fracture - Structure containing information about fracture networks,
%              independent fractures and their conductivity towards the
%              matrix cells they penetrate. See processFracture2D.
%
% OPTIONAL PARAMETERS:
%
%   assemblyType  - Possible values: 1, 2 or 3.
%                   1. Partition fracture grid based on a minimum and
%                      mean desired cell size as specified by the user.
%                   2. Partition fracture grid by specifying a cell size
%                      ratio (<=1) with respect to total length. ex: If
%                      cell size ratio is 0.1, then each fracture will be
%                      partitioned into ~1/0.1=10 parts. Not a good option
%                      when individual fracture length varies by orders of
%                      magnitude.
%                   3. 1 fracture cell for every matrix cell that fracture
%                      penetrates. 
%
%   min_size      - Scalar quantity, useful only if assemblyType = 1 or 2.
%                   For assemblyType = 1, min_size specifies the minimum
%                   cell size the fracture fine grid can have. The function
%                   will override this option if there exist fractures
%                   which are smaller than min_size. For assemblyType = 2,
%                   min_size is used to determine maximum number of
%                   divisions (~1/min_size) for a fracture line.
%
%   cell_size     - Scalar quantity, useful only if assemblyType = 1 or 2.
%                   For assemblyType = 1, cell_size specifies the average
%                   cell size in the fracture fine grid. For assemblyType =
%                   2, cell_size is used to determine average number of
%                   divisions (~1/cell_size) for a fracture line.
%
% RETURNS:
%   F - Structure with the following sub-structures per fracture line:
%       (a) nodes - contains the following subfields
%           -   start - start index w.r.t. global node numbers for the
%               corresponding fracture line
%           -   coords - node coordinates for the corresponding line
%       (b) cells - contains the following subfields
%           -   start - start index w.r.t. global cell numbers for the
%               corresponding fracture line
%           -   num - Number of grid cells in the corresponding fracture
%               line
%
%   fracture - structure with the added structure - 'intersections' which
%              contains the following fields:
%              (a) lines - n-by-2 matrix of pairs of intersecting lines
%              where n is the total number of fracture intersections.
%              (b) coords - coordinates of intersection correspondng to
%              field lines
%
%
% NOTE: 
%   If there is more than 1 fracture 'intersection' within 1 matrix fine
%   cell, do not use assemblyType = 3.
%
% SEE ALSO:
%   processFracture2D, getIndepNetwork, markcells2D

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

F = struct; % Fracture structure containing nodes/faces
node_count = G.nodes.num;
cell_count = G.cells.num;
if numel(fracture.network)<numel(fracture.lines)
    fracture.intersections.lines = [];
    fracture.intersections.coords = [];
end
opt = struct('assemblyType',1,'min_size', 1, 'cell_size', 2);
opt = merge_options(opt, varargin{:});
type = opt.assemblyType;
if type == 1
    assert(isnumeric(opt.min_size) && isnumeric(opt.cell_size),...
        'Specify a numeric value for cell size');
elseif type == 2
    assert(isnumeric(opt.min_size) && isnumeric(opt.cell_size),...
        'Specify a numeric value for partitioning ratio.');
    assert((opt.min_size>0 && opt.min_size<1) && (opt.cell_size>0 && opt.cell_size<1),...
        'Partitioning ratio must be a number between 0 and 1');
end

%-------------------------------------------------------------------------%

if type == 1
    possible_ratios = 0:0.0001:1;
    possible_ratios = possible_ratios(mod(1,possible_ratios)==0);
    for i = 1:numel(fracture.lines)
        min_size = opt.min_size;
        cell_size = opt.cell_size;
        F(i).nodes.start = node_count+1;
        F(i).cells.start = cell_count+1;
        nw = fracture.lines(i).network;
        ol = fracture.network(nw).lines; ol = ol(ol~=i);
        fracp1 = fracture.lines(i).endp; A=fracp1(1:2); B = fracp1(3:4);
        frac_length = pdist_euclid([A;B]);
        ratio = cell_size/frac_length;
        ratio_min = min_size/frac_length;
        [~,loc] = min(abs(possible_ratios-ratio));
        ratio = possible_ratios(loc);
        [~,loc2] = min(abs(possible_ratios-ratio_min));
        if ratio<ratio_min,
            if loc == numel(possible_ratios)
                disp(['Line ',num2str(i),' is too small. It is considered as 1 fracture cell']);
            else
                while(ratio<ratio_min)
                    ratio = possible_ratios(loc+1);
                    if loc+1 >= numel(possible_ratios)
                        disp(['Line ',num2str(i),' is too small. It is considered as 1 fracture cell']);
                    end
                    loc = loc+1;
                end
            end
        end
        cell_size = ratio;
        min_size = possible_ratios(loc2);
        clear loc loc2 ratio ratio_min frac_length
        if ~isempty(ol) % i.e there is > 1 fracture in network
            fracp2 = zeros(numel(ol),G.griddim*2); % 2 end points per fracture line other than primary fracture
            for j = 1:numel(ol)
                fracp2(j,:) = fracture.lines(ol(j)).endp;
            end
            out = lineSegmentIntersect(fracp1,fracp2);
            n_intersect = numel(find(out.intAdjacencyMatrix)); %n_intersect has to be > 1 since 'ol' was not empty
            xi = out.intMatrixX(out.intAdjacencyMatrix).';
            yi = out.intMatrixY(out.intAdjacencyMatrix).';
            fracture.intersections.lines(end+1:end+n_intersect,:) = ...
                sort([repmat(i,n_intersect,1),ol(out.intAdjacencyMatrix).'],2);
            fracture.intersections.coords(end+1:end+n_intersect,:) = [xi,yi];
            tr = zeros(n_intersect,1);
            for j = 1:n_intersect
                tr(j) = abs((xi(j)-B(1))/(A(1)-B(1)));
                if isnan(tr(j)), tr(j) = abs((yi(j)-B(2))/(A(2)-B(2)));end % line // to y-axis
            end
            t = unique([transpose(0:cell_size:1);tr]);
            tdiff = abs(diff(t));
            if any(tdiff-min_size<-1e-15)
                pop = find(tdiff-min_size<-1e-15); % pop these elements such that the fracture intersection is still a node
                ignore = [1;numel(t)];
                while (~isempty(pop))
                    flag = 0;
                    if ismember(t(pop(1)),tr) && ismember(t(pop(1)+1),tr)
                        flag = 1;
                        ignore = [ignore;pop(1)]; %#ok
                    elseif ismember(t(pop(1)),tr)
                        t = [t(1:pop(1));t(pop(1)+2:end)];
                    elseif ismember(t(pop(1)+1),tr)
                        t = [t(1:pop(1)-1);t(pop(1)+1:end)];
                    else error('Missing fracture end points.');
                    end
                    tdiff = abs(diff(t));
                    pop = find(tdiff-min_size<-1e-15);
                    if flag == 1
                        pop = pop(~ismember(pop,ignore));
                    end
                end
            end
            points = [A(1)*t+(1-t)*B(1),A(2)*t+(1-t)*B(2)];
            ncount = size(points,1);
            F(i).nodes.coords = points;
            F(i).cells.num = ncount-1;
        else
            t = transpose(0:cell_size:1);
            points = [A(1)*t+(1-t)*B(1),A(2)*t+(1-t)*B(2)];
            ncount = size(points,1);
            F(i).nodes.coords = points;
            F(i).cells.num = ncount-1;
        end
        node_count = node_count + ncount;
        cell_count = cell_count + ncount - 1;
    end
    if isfield(fracture,'intersections')
        [fracture.intersections.lines,ind] = ...
            unique(fracture.intersections.lines,'rows');
        fracture.intersections.coords = fracture.intersections.coords(ind,:);
    end
    
    %-------------------------------------------------------------------------%
    
elseif type == 2
    min_size = opt.min_size; % 1 must be divisible by this number
    cell_size = opt.cell_size; % 1 must be divisible by this number
    if mod(1,min_size)~=0 || mod(1,cell_size)~=0,error('Select correct cell size ratio'),end
    
    for i = 1:numel(fracture.lines)
        F(i).nodes.start = node_count+1;
        F(i).cells.start = cell_count+1;
        nw = fracture.lines(i).network;
        ol = fracture.network(nw).lines; ol = ol(ol~=i);
        fracp1 = fracture.lines(i).endp; A=fracp1(1:2); B = fracp1(3:4);
        if isempty(ol)==0 % i.e there is > 1 fracture in network
            fracp2 = zeros(numel(ol),G.griddim*2); % 2 end points per fracture line other than primary fracture
            for j = 1:numel(ol)
                fracp2(j,:) = fracture.lines(ol(j)).endp;
            end
            out = lineSegmentIntersect(fracp1,fracp2);
            n_intersect = numel(find(out.intAdjacencyMatrix)); %n_intersect has to be > 1 since 'ol' was not empty
            xi = out.intMatrixX(out.intAdjacencyMatrix).';
            yi = out.intMatrixY(out.intAdjacencyMatrix).';
            fracture.intersections.lines(end+1:end+n_intersect,:) = ...
                sort([repmat(i,n_intersect,1),ol(out.intAdjacencyMatrix).'],2);
            fracture.intersections.coords(end+1:end+n_intersect,:) = [xi,yi];
            tr = zeros(n_intersect,1);
            for j = 1:n_intersect
                tr(j) = abs((xi(j)-B(1))/(A(1)-B(1)));
                if isnan(tr(j)), tr(j) = abs((yi(j)-B(2))/(A(2)-B(2)));end % line // to y-axis
            end
            t = unique([transpose(0:cell_size:1);tr]);
            tdiff = abs(diff(t));
            if any(tdiff-min_size<-1e-15)
                pop = find(tdiff-min_size<-1e-15); % pop these elements such that the fracture intersection is still a node
                while isempty(pop)==0
                    if ismember(t(pop(1)),tr)
                        t = [t(1:pop(1));t(pop(1)+2:end)];
                    elseif ismember(t(pop(1)+1),tr)
                        t = [t(1:pop(1)-1);t(pop(1)+1:end)];
                    else error('Missing fracture end points.');
                    end
                    tdiff = abs(diff(t));
                    pop = find(tdiff-min_size<-1e-15);
                end
            end
            points = [A(1)*t+(1-t)*B(1),A(2)*t+(1-t)*B(2)];
            ncount = size(points,1);
            F(i).nodes.coords = points;
            F(i).cells.num = ncount-1;
        else
            t = transpose(0:cell_size:1);
            points = [A(1)*t+(1-t)*B(1),A(2)*t+(1-t)*B(2)];
            ncount = size(points,1);
            F(i).nodes.coords = points;
            F(i).cells.num = ncount-1;
        end
        node_count = node_count + ncount;
        cell_count = cell_count + ncount - 1;
    end
    if isfield(fracture,'intersections')
        [fracture.intersections.lines,ind] = ...
            unique(fracture.intersections.lines,'rows');
        fracture.intersections.coords = fracture.intersections.coords(ind,:);
    end

%-------------------------------------------------------------------------%
  
else
    for i = 1:numel(fracture.network) % Nodes and cells defined per network.
        for j = 1:numel(fracture.network(i).lines)
            fl = fracture.network(i).lines(j);
            fcells = fracture.lines(fl).cells; % Fine cells belonging to the particular fracture line
            frac_endp = fracture.lines(fl).endp;
            ccount = 0;
            ncount = 0;
            F(fl).nodes.start = node_count+1;
            F(fl).cells.start = cell_count+1;
            coords = [];
            for k = 1:numel(fcells)
                ci = fcells(k); % Fine cell containing fracture
                
                % faces and nodes for cell - ci
                faces = G.cells.faces(G.cells.facePos(ci):G.cells.facePos(ci+1)-1,1);
                edges = zeros(numel(faces),4);
                cnodes = [];
                for l = 1:numel(faces)
                    fnodes = G.faces.nodes(G.faces.nodePos(faces(l)):G.faces.nodePos(faces(l)+1)-1); % 2 for 2D, 4 for 3D
                    cnodes = [cnodes;fnodes];%#ok
                    edges(l,:) = reshape(G.nodes.coords(fnodes,:).',1,numel(G.nodes.coords(fnodes,:))); % faces represented as line segments
                end
                cnodes = unique(cnodes);
                out = lineSegmentIntersect(frac_endp,edges);
                xi = out.intMatrixX(out.intAdjacencyMatrix).';
                yi = out.intMatrixY(out.intAdjacencyMatrix).';
                
                % Now, check for the extent of penetration of fracture and
                % intersection with any other fractures inside the cell
                flag = 0; % Only 1 fracture goes through the cell intersecting at 2 faces
                
                if numel(find(out.intAdjacencyMatrix))==1
                    flag = 1;
                elseif numel(find(out.intAdjacencyMatrix))==2
                    if diff(xi)+diff(yi)==0
                        flag = 1;
                    end
                end
                % flag = 1 if only 1 fracture intersects at 1 face only.
                
                if numel([G.cells.fracture.line_num{ci,1},G.cells.fracture.dfm_line_num{ci,1}])>1 % More than 1 fracture in the cell
                    % Check if they intersect first. Option 1 does not yet deal
                    % with multiple intersections in 1 cell. If that is the
                    % case use Option 2
                    fl2 = [G.cells.fracture.line_num{ci,1},G.cells.fracture.dfm_line_num{ci,1}];
                    fl2 = fl2(fl2~=fl);
                    fracp2 = fracture.lines(fl2).endp; % Can only be 1 fracture
                    if flag == 1
                        points = [G.cells.centroids(ci,:);G.nodes.coords(cnodes,:)];
                        tri = delaunayTriangulation(points(:,1),points(:,2));
                        K = convexHull(tri);
                        fendp = [frac_endp(1:2);frac_endp(3:4)];
                        [in,~] = inpolygon(fendp(:,1),fendp(:,2),tri.Points(K,1),tri.Points(K,2));
                        if find(in)==1 % fracture starting point
                            fracp = [fendp(in,:); uniquetol(xi,eps), uniquetol(yi,eps)];
                        else % fracture end point
                            fracp = [uniquetol(xi,eps), uniquetol(yi,eps);fendp(in,:)];
                        end
                        out = lineSegmentIntersect([fracp(1,:),fracp(2,:)],fracp2);
                        if out.intAdjacencyMatrix % Intersection=true
                            flag = 3;
                        end
                    else
                        fracp = uniquetol([xi,yi],eps,'ByRows',true);
                        out = lineSegmentIntersect([fracp(1,:),fracp(2,:)],fracp2);
                        if out.intAdjacencyMatrix % Intersection = true
                            flag = 2;
                        end
                    end
                end
                %%
                switch flag
                    case 0 % 1 fracture in cell cutting at 2 faces (i.e. going through)
                        fracp = uniquetol([xi,yi],eps,'ByRows',true);
                        ccount = ccount+1;
                        coords(ncount+1:ncount+2,:) = fracp;
                        ncount = ncount+2;
                    case 1 % 1 fracture in cell only partially penetrating (i.e. fracture end-point in cell)
                        points = [G.cells.centroids(ci,:);G.nodes.coords(cnodes,:)];
                        tri = delaunayTriangulation(points(:,1),points(:,2));
                        K = convexHull(tri);
                        fendp = [frac_endp(1:2);frac_endp(3:4)];
                        [in,~] = inpolygon(fendp(:,1),fendp(:,2),tri.Points(K,1),tri.Points(K,2));
                        fracp = [fendp(in,:); uniquetol(xi,eps), uniquetol(yi,eps)];
                        ccount = ccount+1;
                        coords(ncount+1:ncount+2,:) = fracp;
                        ncount = ncount+2;
                    case 2 % 2 fractures intersecting in cell
                        xi = out.intMatrixX(out.intAdjacencyMatrix).';
                        yi = out.intMatrixY(out.intAdjacencyMatrix).';
                        fracture.intersections.lines(end+1,:) = unique([fl,fl2]);
                        fracture.intersections.coords(end+1,:) = [xi,yi];
                        coords(ncount+1:ncount+2,:) = [fracp(1,:);xi,yi];
                        coords(ncount+3:ncount+4,:) = [xi,yi;fracp(2,:)];
                        ccount = ccount+2;
                        ncount = ncount+4;
                    case 3 % 2 fractures intersecting with 1 or both partially penetrating (i.e. end points in cell)
                        xi = out.intMatrixX(out.intAdjacencyMatrix).';
                        yi = out.intMatrixY(out.intAdjacencyMatrix).';
                        fracture.intersections.lines(end+1,:) = unique([fl,fl2]);
                        fracture.intersections.coords(end+1,:) = [xi,yi];
                        coords(ncount+1:ncount+2,:) = [fracp(1,:);xi,yi];
                        if isequal([xi,yi],fracp(2,:))==0
                            coords(ncount+3:ncount+4,:) = [xi,yi;fracp(2,:)];
                            ncount = ncount+4;
                            ccount = ccount+2;
                        else % fracture line end point is the intersection itself
                            ncount = ncount+2;
                            ccount = ccount+1;
                        end
                end
            end
            coords = uniquetol(coords,eps,'ByRows',true);
            F(fl).nodes.coords = coords;
            ncells = size(coords,1)-1;
            if ncells~=ccount && ncells~=ccount/2, warning('Possible error in indexing'), end
            F(fl).cells.num = ncells;
            cell_count = cell_count+ncells;
            node_count = node_count+ncells+1;
        end
        if isfield(fracture,'intersections')
            [fracture.intersections.lines,ind] = ...
                unique(fracture.intersections.lines,'rows');
            fracture.intersections.coords = fracture.intersections.coords(ind,:);
        end
    end
end

%% Rearrange frac nodes to pass monotonically varying coordinates into tensorGrid

for i = 1:numel(F)
    coords = F(i).nodes.coords;
    % uniquetol(coords,eps,'ByRows',true);
    endp = [coords(1,:);coords(end,:)];
    diff1 = diff(endp,1);
    if abs(diff1(2)) < eps*100 % // to x-axis
        [~,ind] = sort(coords(:,1)); % Sort starting from lowest x-coord
        coords = coords(ind,:);
    elseif abs(diff1(1)) < eps*100 % // to y-axis
        [~,ind] = sort(coords(:,2)); % Sort starting from lowest x-coord
        coords = coords(ind,:);
    else
        [~,ind] = sort(coords(:,1)); % Sort starting from lowest x-coord
        coords = coords(ind,:);
    end
    F(i).nodes.coords = coords;
end

return