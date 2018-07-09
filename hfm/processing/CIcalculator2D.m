function G = CIcalculator2D(G,fracture,varargin)
% CIcalculator2D computes the conductivity index (CI) of each 2D cell for
% every fracture line embedded in it.
%
% SYNOPSIS:
%   G = CIcalculator2D(G,fracture)
%   G = CIcalculator2D(G,fracture, waitbar)
%
% DESCRIPTION:
%   CI or conductivity index is similar to the well productivity index as
%   in the Peaceman well model. In hierarchical/embedded fracture
%   modelling, the flux exchange between fracture and matrix is defined by
%   a model similar to the peaceman well model. See Lee et al, Water
%   Resources Research, 2001 or SPE-65095-PA, Lee et al, 2000.
%
% REQUIRED PARAMETERS:
%
%   G         - Matrix grid structure as returned by processFracture2D.
%
%   fracture  - Processed fracture structure. See processFracture2D.
%
%   waitbar   - Display a waitbar if this boolean variable is true
%
% RETURNS:
%   G  - Matrix grid structure with added cell lists 'G.cells.fracture.CI'
%        and 'G.cells.fracture.fA' containing fracture-matrix conductivity
%        index and fracture length inside each matrix cell respectively
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

if nargin>2
   showProgress=varargin{1};
   hwb=waitbar(0,'Computing conductivity index (embedded)');
else
   showProgress=false;
end
frac_cells = find(G.cells.fracture.indicator==1); % Matrix cells containing embedded fractures
CI = cell(G.cells.num,1);
fA = cell(G.cells.num,1);
for i = 1:numel(frac_cells)
    flines = cell2mat(G.cells.fracture.line_num(frac_cells(i))); % Fracture lines
    for j = 1:numel(flines)
        fli = flines(j);
        frac_endp = fracture.lines(fli).endp;
        endp_dst = pdist_euclid([frac_endp(1:2);frac_endp(3:4)]);
        faces = G.cells.faces(G.cells.facePos(frac_cells(i)):G.cells.facePos(frac_cells(i)+1)-1,1);
        edges = zeros(numel(faces),4);
        cnodes = [];
        for k = 1:numel(faces)
            fnodes = G.faces.nodes(G.faces.nodePos(faces(k)):G.faces.nodePos(faces(k)+1)-1); % 2 for 2D, 4 for 3D
            cnodes = [cnodes;fnodes]; %#ok
            edges(k,:) = reshape(G.nodes.coords(fnodes,:).',1,numel(G.nodes.coords(fnodes,:))); % faces represented as line segments
        end
        cnodes = unique(cnodes);
        out = lineSegmentIntersect(frac_endp,edges);
        xi = out.intMatrixX(out.intAdjacencyMatrix).';
        yi = out.intMatrixY(out.intAdjacencyMatrix).';
        nc = G.nodes.coords(cnodes,:);
        for failcheck = 1:numel(cnodes)
            d1 = pdist_euclid([frac_endp(1:2);nc(failcheck,:)]);
            d2 = pdist_euclid([frac_endp(3:4);nc(failcheck,:)]);
            if abs((d1+d2)-endp_dst)<eps
               if ~ismember(roundsd(nc(failcheck,:),14),roundsd([xi,yi],14),'rows')
                  xi = [xi;nc(failcheck,1)]; %#ok
                  yi = [yi;nc(failcheck,2)]; %#ok
               end
            end
        end
        points = [xi,yi];
        [~,ia] = unique(roundsd(points,14),'rows');
        points = points(ia,:);
        xi = points(:,1); yi = points(:,2);
        flag = 0;
        if size(points,1)==1
            flag = 1;
        elseif size(points,1)==2
            if diff(xi)+diff(yi)==0 && diff(xi)-diff(yi)==0
                flag = 1;
            end
        end
        if flag == 1
            [new_endp,ratio] = extend_frac2D(G, frac_endp,edges,out,frac_cells(i),cnodes,faces);
            fracp = [ [xi, yi]; new_endp];
            d_avg = getAvgFracDist2D(G, fracp, frac_cells(i), cnodes);
        else
            fracp = [xi,yi];
            [~,ia] = unique(roundsd(fracp,14),'rows');
            fracp = fracp(ia,:);
            d_avg = getAvgFracDist2D(G, fracp, frac_cells(i), cnodes);
            ratio = 1;
        end
        CI{frac_cells(i),1} = [CI{frac_cells(i),1}, (pdist_euclid(fracp)/d_avg)*ratio];
        fA{frac_cells(i),1} = [fA{frac_cells(i),1}, pdist_euclid(fracp)*ratio];
    end
    if showProgress,
       waitbar(i/numel(frac_cells),hwb);
    end
end
if showProgress
   close(hwb);
   hwb=waitbar(0,'Computing conductivity index (discrete)');
end
frac_cells = find(G.cells.fracture.dfm_indicator==1);
for i = 1:numel(frac_cells)
    flines = cell2mat(G.cells.fracture.dfm_line_num(frac_cells(i)));
    for j = 1:numel(flines)
        fli = flines(j);
        frac_endp = fracture.lines(fli).endp;
        endp_dst = pdist_euclid([frac_endp(1:2);frac_endp(3:4)]);
        faces = G.cells.faces(G.cells.facePos(frac_cells(i)):G.cells.facePos(frac_cells(i)+1)-1,1);
        edges = zeros(numel(faces),4);
        cnodes = [];
        for k = 1:numel(faces)
            fnodes = G.faces.nodes(G.faces.nodePos(faces(k)):G.faces.nodePos(faces(k)+1)-1);
            cnodes = [cnodes;fnodes]; %#ok
            edges(k,:) = reshape(G.nodes.coords(fnodes,:).',1,numel(G.nodes.coords(fnodes,:)));
        end
        cnodes = unique(cnodes);
        out = lineSegmentIntersect(frac_endp,edges);
        xi = out.intMatrixX(out.intAdjacencyMatrix).';
        yi = out.intMatrixY(out.intAdjacencyMatrix).';
        nc = G.nodes.coords(cnodes,:);
        for failcheck = 1:numel(cnodes)
            d1 = pdist_euclid([frac_endp(1:2);nc(failcheck,:)]);
            d2 = pdist_euclid([frac_endp(3:4);nc(failcheck,:)]);
            if abs((d1+d2)-endp_dst)<eps
               if ~ismember(roundsd(nc(failcheck,:),14), roundsd([xi,yi],14), 'rows')
                  xi = [xi;nc(failcheck,1)]; %#ok
                  yi = [yi;nc(failcheck,2)]; %#ok
               end
            end
        end
        points = [xi,yi];
        [~,ia] = unique(roundsd(points,14),'rows');
        points = points(ia,:);
        xi = points(:,1); yi = points(:,2);
        flag = 0;
        if size(points,1)==1
            flag = 1;
        elseif size(points,1)==2
            if diff(xi)+diff(yi)==0 && diff(xi)-diff(yi)==0
                flag = 1;
            end
        end
        if flag==1
            xi = unique(out.intMatrixX(out.intAdjacencyMatrix));
            yi = unique(out.intMatrixY(out.intAdjacencyMatrix));
            coinc_edge = edges(out.coincAdjacencyMatrix,:);
%             if isempty(coinc_edge) % temp fix - reduce significant digits
%                 temp = find(abs(out.intNormalizedDistance2To1)<eps*100);
%                 coinc_edge = edges(temp(1)+1,:);
%             end
            dist = pdist_euclid([coinc_edge(1:2);coinc_edge(3:4)]);
            dist1 = pdist_euclid([frac_endp(1:2);xi,yi]);
            dist2 = pdist_euclid([frac_endp(3:4);xi,yi]);
            if dist1<=dist2
                ratio = dist1/dist;
            else
                ratio = dist2/dist;
            end
            % Rearrange from [x1 y1 x2 y2] to [x1 y1; x2 y2] format
            fracp = [coinc_edge(1:2);coinc_edge(3:4)];
            d_avg = getAvgFracDist2D(G, fracp, frac_cells(i), cnodes);
        else
            coinc_edge = edges(out.coincAdjacencyMatrix,:);
%             if isempty(coinc_edge) % temp fix - reduce significant digits
%                 temp = find(abs(out.intNormalizedDistance2To1)<eps*100);
%                 coinc_edge = edges(temp(1)+1,:);
%             end
            fracp = [coinc_edge(1:2);coinc_edge(3:4)];
            d_avg = getAvgFracDist2D(G, fracp, frac_cells(i), cnodes);
            ratio = 1;
        end
        CI{frac_cells(i),1} = [CI{frac_cells(i),1}, (pdist_euclid(fracp)/d_avg)*ratio];
        fA{frac_cells(i),1} = [fA{frac_cells(i),1}, pdist_euclid(fracp)*ratio];
    end
    if showProgress,
       waitbar(i/numel(frac_cells),hwb);
    end
end
if showProgress,
   close(hwb);
end
G.cells.fracture.CI = CI;
G.cells.fracture.fA = fA;
return
        