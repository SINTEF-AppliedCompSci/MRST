function G = frac_matrix_nnc(G,F,fracture)
% frac_matrix_nnc assigns a "non-neighboring connection (NNC)" indicators
% to each fracture-matrix connection and also assigns a transmissibility to
% each NNC. See Lee et al, Water Resources Research, 2001 or SPE-65095-PA,
% Lee et al, 2000.
%
% SYNOPSIS:
%   G = frac_matrix_nnc(G,F,fracture)
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing G.FracGrid (see
%                 FracTensorGrid2D) and G.cells.fracture (see markcells2D
%                 and CIcalculator2D)
%
%   F, fracture - Output from gridFracture2D.
%
% RETURNS:
%   G - Grid structure with fracture-matrix connections  at fine scale and
%       their corresponding CI added in 'G.nnc' as lists 'G.nnc.cells' and
%       'G.nnc.CI' respectively. To aid in detecting specific NNC types,
%       these connections are added as 'frac-matrix' type in the list
%       G.nnc.type
%
% SEE ALSO:
%   getIndepNetwork, markcells2D, CIcalculator2D, gridFracture2D,
%   frac_frac_nnc 

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


G.nnc.cells = zeros(0,2);
G.nnc.CI = zeros(0,1);
count = 1;
for i = 1:numel(F)
    mcells = fracture.lines(i).cells; % Matrix cells which contain fracture line 1
    for j = 1:numel(mcells)
        clear node
        CI_tot = G.cells.fracture.CI{mcells(j),1};
        lines_in_cell = [G.cells.fracture.line_num{mcells(j),1},...
            G.cells.fracture.dfm_line_num{mcells(j),1}];
        CI_tot = CI_tot(lines_in_cell==i);
        fA = G.cells.fracture.fA{mcells(j),1};
        fA = fA(lines_in_cell==i); % fracture length in cell
        % Assemble cell nodes and edges
        faces = G.cells.faces(G.cells.facePos(mcells(j)):G.cells.facePos(mcells(j)+1)-1,1);
        edges = zeros(numel(faces),4);
        cnodes = [];
        for k = 1:numel(faces)
            fnodes = G.faces.nodes(G.faces.nodePos(faces(k)):G.faces.nodePos(faces(k)+1)-1); % 2 for 2D, 4 for 3D
            cnodes = [cnodes;fnodes]; %#ok
            edges(k,:) = reshape(G.nodes.coords(fnodes,:).',1,numel(G.nodes.coords(fnodes,:))); % faces represented as line segments
        end
        cnodes = unique(cnodes);
        points = [G.cells.centroids(mcells(j),:);G.nodes.coords(cnodes,:)];
        tri = delaunayTriangulation(points(:,1),points(:,2));
        K = convexHull(tri);
        xq = F(i).nodes.coords(:,1);
        yq = F(i).nodes.coords(:,2);
        [in, on] = inpolygon(xq,yq,tri.Points(K,1),tri.Points(K,2));
        
        if (any(in)==0 && any(on)) || (any(in) && isequal(in,on))
            % In this case 1 matrix cell is interacting with 1 fracture
            % cell.
            %-------------------------------------------------------------%
            % No nodes inside the cell, they could be on the face of the
            % cell or the fracture nodes may not be inside this cell
            % altogether. This could happen when a fracture grid cell is
            % much larger than the section of the fracture that penetrates
            % this particular cell.
            
            if numel(find(on))==1
                % When there is just one node on the cell edge.
                node = find(on);
                if node == 1
                    fcellno = F(i).cells.start + node - 1;
                elseif node == size(F(i).nodes.coords,1)
                    fcellno = F(i).cells.start + node - 2;
                else
                    nprev = F(i).nodes.coords(node-1:node,:);
                    nprev = [nprev(1,:), nprev(2,:)];
                    nnext = F(i).nodes.coords(node:node+1,:);
                    nnext = [nnext(1,:), nnext(2,:)];
                    oprev = lineSegmentIntersect(nprev, edges);
                    onext = lineSegmentIntersect(nnext, edges);
                    pprev = unique([oprev.intMatrixX(oprev.intAdjacencyMatrix).',...
                        oprev.intMatrixY(oprev.intAdjacencyMatrix).'],'rows');
                    pnext = unique([onext.intMatrixX(onext.intAdjacencyMatrix).',...
                        onext.intMatrixY(onext.intAdjacencyMatrix).'],'rows');
                    if size(pprev,1)>size(pnext,1)
                        fcellno = F(i).cells.start + node - 2;
                    elseif size(pprev,1)<size(pnext,1)
                        fcellno = F(i).cells.start + node - 1;
                    else
                        fcellno(1) = F(i).cells.start + node - 2;
                        fcellno(2) = F(i).cells.start + node - 1;
                        tdist = pdist_euclid([pprev;pnext]);
                        fpts = [xq(on),yq(on)];
                        r(1) = pdist_euclid([pprev;fpts])/tdist;
                        r(2) = pdist_euclid([pnext;fpts])/tdist;
                        for k = 1:2
                            G.nnc.cells(count,:) = [mcells(j) fcellno(k)];
                            G.nnc.CI(count,1) = CI_tot*r(k);
                            count = count+1;
                        end
                    end
                end
                if numel(fcellno)==1
                    G.nnc.cells(count,:) = [mcells(j) fcellno];
                    G.nnc.CI(count,1) = CI_tot;
                    count = count+1;
                end 
            elseif numel(find(on))==2
                % When the fracture segment inside the cell is a fracture grid
                % cell itself
                node = find(on);
                fcellno = F(i).cells.start + node(1) - 1;
                G.nnc.cells(count,:) = [mcells(j) fcellno];
                G.nnc.CI(count,1) = CI_tot;
                count = count+1;
            end
        elseif any(in+on)==0
            % Fracture grid size is much larger than the length of the
            % segment inside this matrix cell
            for k = 1:size(F(i).nodes.coords,1)-1
                ncoords = [F(i).nodes.coords(k,:), F(i).nodes.coords(k+1,:)];
                out = lineSegmentIntersect(ncoords, edges);
                if any(out.intAdjacencyMatrix)
                    G.nnc.cells(count,:) = [mcells(j) k+F(i).cells.start-1];
                    G.nnc.CI(count,1) = CI_tot;
                    count = count+1;
                end
            end
        else
            node = find(in&~on); % only ones inside cell
            flag = 0;
            %-----------------------------REDO----------------------------%
            if node(1) ~= numel(in) && node(1)~=1 && ~any(node == numel(in))  % Fracture passes through this cell
                flag = 1;
                node_new = [node(1)-1;node]; % fracture grid cell numbers 
                node = node_new; clear node_new
            elseif any(node == numel(in)) && numel(node)>1 % Fracture ends "inside" this cell with >1 nodes in cell
                flag = 3;
                node_new = [node(1)-1;node];
                node = node_new; clear node_new 
            elseif node(1) == 1 && numel(node)>1 % Fracture starts "inside" this cell with >1 nodes in cell
                flag = 2;
                node_new = [node;node(end)+1];
                node = node_new; clear node_new
            elseif node == numel(in) % Only fracture endp in cell
                node = node-1;
            end
            node = unique(node); % because of 2nd elseif
            %----Make the NNC list --> [Matrix cell ind, Frac cell ind]---%
            if flag == 3
                G.nnc.cells(count:count+numel(node)-2,:) = ...
                    [repmat(mcells(j),numel(node)-1,1),node(1:end-1)+F(i).cells.start-1];
            elseif flag == 2
                G.nnc.cells(count:count+numel(node)-2,:) = [repmat(mcells(j),...
                    numel(node)-1,1),node(1:end-1)+F(i).cells.start-1];
            else
                G.nnc.cells(count:count+numel(node)-1,:) = [repmat(mcells(j),...
                    numel(node),1),node+F(i).cells.start-1];
            end
            if numel(node)==1 
                % May happen at fracture end points or could also happen if
                % fracture resolution is very high but in that case the
                % previous statements modify node to ensure numel(node)~=1.
                G.nnc.CI(count,1) = CI_tot;
            else
                ncoords = [F(i).nodes.coords(node(1:end-1),:), F(i).nodes.coords(node(2:end),:)];
                if flag == 1 || flag == 3
                    out = lineSegmentIntersect(ncoords(1,:), edges);
                    pi = unique([out.intMatrixX(out.intAdjacencyMatrix).',...
                        out.intMatrixY(out.intAdjacencyMatrix).'],'rows');
                    if isempty(pi)
                        if any(abs(out.intNormalizedDistance1To2)<1e-10)
                        % For the unlikely case that lineSegmentIntersect
                        % does not detect a point on the edges 
                            points = F(i).nodes.coords(node,:);
                        else
                            ind = closestNode(ncoords,G.nodes.coords(cnodes,:));
                            points = [G.nodes.coords(cnodes(ind),:);F(i).nodes.coords(node(2:end),:)];
                        end
                    else
                        points = [pi;F(i).nodes.coords(node(2:end),:)];
                    end
                elseif flag == 2
                    out = lineSegmentIntersect(ncoords(end,:), edges);
                    pi = unique([out.intMatrixX(out.intAdjacencyMatrix).',...
                        out.intMatrixY(out.intAdjacencyMatrix).'],'rows');
                    if isempty(pi) 
                        if any(abs(out.intNormalizedDistance2To1)<1e-10)
                        % For the unlikely case that lineSegmentIntersect
                        % does not detect a point on the edges
                            points = F(i).nodes.coords(node,:);
                        else
                            ind = closestNode(ncoords,G.nodes.coords(cnodes,:));
                            points = [F(i).nodes.coords(node(2:end),:);G.nodes.coords(cnodes(ind),:)];
                        end    
                    else
                        points = [F(i).nodes.coords(node(1:end-1),:);pi];
                    end
                end
                points = unique(points,'rows','stable');
                if flag == 3
                    ratio = zeros(numel(node)-1,1);
                    for k = 2:numel(node)-1
                        ratio(k-1) = pdist_euclid(points(k-1:k,:))/fA;
                    end
                elseif flag ==2
                    ratio = zeros(numel(node)-1,1);
                    for k = 1:numel(node)-2
                        ratio(k) = pdist_euclid(points(k:k+1,:))/fA;
                    end
                else
                    ratio = zeros(numel(node),1);
                    for k = 2:numel(node)
                        ratio(k-1) = pdist_euclid(points(k-1:k,:))/fA;
                    end
                end
                ratio(end) = 1-sum(ratio);
                if any(ratio<0)
                    if abs(ratio(ratio<0)-0)>1e-6
                        error('Fracture nodes not recognized, negative NNC strength detected.');
                    else
                        % Ignore very small fracture-matrix connections.
                        ratio(ratio<0) = 0;
                    end
                end
                if flag == 3 || flag == 2
                    G.nnc.CI(count:count+numel(node)-2,1) = ratio*CI_tot;
                else
                    G.nnc.CI(count:count+numel(node)-1,1) = ratio*CI_tot;
                end
            end
            count = count+numel(node);
        end
        clear faces edges points
    end
end
G.nnc.type = repmat({'frac-matrix'},numel(G.nnc.CI),1);
return

function ni = closestNode(endp,npoints)
%
endp = [endp(1:2);endp(3:4)];
slope = diff(endp(:,2))/diff(endp(:,1));
if isinf(slope)
    A = 1; B = 0; C = -endp(1,1);
elseif slope == 0
    A = 0; B = 1; C = -endp(1,2);
else
    A = -slope;
    B = 1;
    C = -endp(1,2)+endp(1,1)*slope;
end
%
d = abs(A*npoints(:,1)+B*npoints(:,2)+C)./(sqrt(A^2+B^2));
[~,mind] = min(d);
if numel(mind) == 1
    ni = mind; return
else
    for nps = 1:numel(mind)
        dsum = pdist_euclid([npoints(mind(nps),:);endp]);
        if abs(sum(dsum(1:2))-dsum(3))<1e-10
            ni = mind(nps);
            return
        end
    end
end
