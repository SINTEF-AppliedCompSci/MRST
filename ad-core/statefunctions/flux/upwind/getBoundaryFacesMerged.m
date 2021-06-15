function [bfaces, bcells, bfc] = getBoundaryFacesMerged(G, dim)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

    if nargin == 1
        dim = G.griddim;
    end
    
    [bfaces, bcells] = boundaryFaces(G);
    bfc = G.faces.centroids(bfaces, 1:dim);
    bcc = G.cells.centroids(bcells, 1:dim);
    nbf = numel(bfaces);
%     % Extrude to mimick ghost cells
    bfc = bcc + 2*(bfc - bcc);
    if 0%any(strcmpi(G.type, 'tensorGrid') | strcmpi(G.type, 'processGRDECL'));
        
        isbf = false(G.faces.num, 1);
        isbf(bfaces) = true;
        facetags = zeros(nbf, 1);
        gcf = G.cells.faces;
        
        for i = 1:nbf
            facetags(i) = gcf(gcf(:, 1) == bfaces(i), 2);
        end
        
        keep = true(nbf, 1);
        for i = 1:nbf
            c = bcells(i);
            
            fpos = G.cells.facePos(c):(G.cells.facePos(c+1)-1);
            
            faces = G.cells.faces(fpos, 1);
            tags = G.cells.faces(fpos, 2);
            
            sub = faces == bfaces(i);
            
            tag = tags(sub);
            
            subf = faces(tags == tag);
            isbfsub = isbf(subf);
%             keep(i) = all(isbfsub);
            keep(i) = all(isbfsub) && subf(1) == bfaces(i);
        end
        
        if dim == 2
            keep = keep & facetags < 5;
        end
        
        bfaces = bfaces(keep);
        bcells = bcells(keep);
        bfc = bfc(keep, :);
        return
    end
    
    return
    % Process to merge overlapping boundaries.
    G = createAugmentedGrid(G);
    
    gcf = G.cells.faces;
    for i = 1:nbf
        tags(i) = gcf(gcf(:, 1) == bfaces(i), 2);
    end
    
    bf_partition = (1:numel(bfaces))';
    
    for f = 1:numel(bcells)
        c = bcells(f);
        
        subf = bcells == c;
        
        if nnz(subf) > 1
            
            if 1
                bfaceNo = find(subf);
%                 faces = bfaces(subf);
                
                ftags = tags(bfaceNo);

                for cno = 1:max(ftags)
                    ff = bfaceNo(ftags == cno);
                    if ~isempty(ff)
                        bf_partition(ff) = bf_partition(ff(1));
                    end
                end         
            else
                if 0
                    ff = find(subf);
                    bf_partition(ff) = bf_partition(ff(1));
                else
                % Find shared edges for boundary faces for this cell
                bfaceNo = find(subf);
                faces = bfaces(subf);

                if 0
                    positions = G.faces.edgePos;
                    values = G.faces.edges;
                    maxCounts = 1;
                else
                    positions = G.faces.nodePos;
                    values = G.faces.nodes;
                    maxCounts = 1;
                end

                ep_start = positions(faces);
                ep_end = positions(faces+1)-1;
                counts = ep_end - ep_start + 1;

                faceNo = rldecode((1:numel(faces))', counts);
                edgePos = mcolon(ep_start, ep_end);
                edges = values(edgePos);


                occurances = accumarray(edges, 1);
                occurances = occurances(edges);



                overlap = occurances > maxCounts;

                if any(overlap)

                    mat = sparse(edges(overlap), faceNo(overlap), 1, max(edges), numel(faces));
                    mat = mat(edges, :);

                    A = mat'*mat;
    %                 A = A + A';
                    cc = components(A);

                    for cno = 1:max(cc)
                        ff = bfaceNo(cc == cno);
                        bf_partition(ff) = bf_partition(ff(1));
                    end
                end
    %             for e = 1:numel(occurances)
    %                 if occurances(e) > 1
    %                     ff = bfaceNo(faceNo(edges == edges(e)));
    %                     bf_partition(ff) = bf_partition(ff(1));
    %                 end
    %             end
                end
            end
        end
    end
    keep = unique(bf_partition);
    if 1
        % Weight each face by area and define new centroids as the weighted
        % sum
        wval = G.faces.areas(bfaces);
%         wval = 1;
        beta = accumarray((1:numel(bfaces))', wval);
        beta_tot = accumarray(bf_partition, wval);

        w = beta./beta_tot(bf_partition);
        
        dim = size(bfc, 2);
        bfc_new = zeros(numel(keep), dim);
        
        for i = 1:dim
            tmp = accumarray(bf_partition, w.*bfc(:, i));
            bfc_new(:, i) = tmp(keep);
        end
        bfc = bfc_new;
    else
        bfc = bfc(keep, :);
    end
    bfaces = bfaces(keep);
    bcells = bcells(keep); 
end
