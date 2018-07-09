function G = frac_frac_nnc(G,F,fracture)
% frac_frac_nnc assigns NNC connections to fracture-fracture intersections
% and also assigns a transmissibility to each connection using the
% star-delta transformation. See see SPE-88812-PA, Karimi-Fard et al, 2004.
%
% SYNOPSIS:
%   G = frac_frac_nnc(G,F,fracture)
%
% REQUIRED PARAMETERS:
%
%   G           - Grid data structure containing G.FracGrid (see
%                 FracTensorGrid2D) and corresponding rock properties for
%                 the fracture.
%
%   F, fracture - Output from gridFracture2D.
%
% RETURNS:
%   G - Grid structure with fracture-fracture intersections and their
%       transmissibilities added in the cell lists 'G.nnc.cells' and
%       'G.nnc.T' respectively. To aid in detecting specific NNC types,
%       these connections are added as 'star-delta' type in the list
%       'G.nnc.type'
%
% SEE ALSO:
%   assembleFracNodes2D, gridFracture2D, frac_matrix_nnc

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

if isfield(fracture,'intersections')
for i = 1:size(fracture.intersections.lines,1) %% 2D for now
    lines = fracture.intersections.lines(i,:);
    coords = fracture.intersections.coords(i,:);
    [~,Gface(1)] = ismember(roundsd(coords,5),roundsd(F(lines(1)).nodes.coords,5),'rows');
    [~,Gface(2)] = ismember(roundsd(coords,5),roundsd(F(lines(2)).nodes.coords,5),'rows');
    diff1 = diff([F(lines(1)).nodes.coords(1,:);F(lines(1)).nodes.coords(end,:)]);
    diff2 = diff([F(lines(2)).nodes.coords(1,:);F(lines(2)).nodes.coords(end,:)]);
    cells_l = cell(2,1);
    for j = 1:numel(Gface)
        if Gface(j) == 1
            cells_l{j,1} = Gface(j);
        elseif Gface(j) == size(F(lines(j)).nodes.coords,1)
            cells_l{j,1} = Gface(j)-1;
        else
            cells_l{j,1} = [Gface(j)-1,Gface(j)];
        end
    end
    Gf1 = G.FracGrid.(['Frac',num2str(lines(1))]);
    Gf2 = G.FracGrid.(['Frac',num2str(lines(2))]);
    cells_l{1,1} = Gf1.cells.start-1 + cells_l{1,1};
    cells_l{2,1} = Gf2.cells.start-1 + cells_l{2,1};

    %-----------Get half face transmissibilities--------------------------%

    T1 = computeTrans(Gf1,Gf1.rock);
    T2 = computeTrans(Gf2,Gf2.rock);
    if diff1(1) == 0
        Gface(1) = Gf1.faces.num - (size(F(lines(1)).nodes.coords,1)-Gface(1));
    end
    if diff2(1) == 0
        Gface(2) = Gf2.faces.num - (size(F(lines(2)).nodes.coords,1)-Gface(2));
    end
    T1 = T1(Gf1.cells.faces(:,1)==Gface(1)); % Half-trans multiplers. t1(1) => for cells_l1(1)
    T2 = T2(Gf2.cells.faces(:,1)==Gface(2)); % Half-trans multiplers. t2(1) => for cells_l2(1)

    %-----------Define NNC's and respective trans using star-delta--------%

    [c1,c2] = meshgrid(cells_l{1,1}, cells_l{2,1});
    [t1,t2] = meshgrid(T1, T2);
    pairs = [c1(:) c2(:)];
    tpairs = [t1(:) t2(:)]; clear c1 c2 t1 t2
    if numel(cells_l{2,1})>1 
        pairs = [cells_l{2,1};pairs]; %#ok
        tpairs = [T2';tpairs]; %#ok
    end
    if numel(cells_l{1,1})>1
        pairs = [cells_l{1,1};pairs]; %#ok
        tpairs = [T1';tpairs]; %#ok
    end
    G.nnc.cells(end+1:end+size(pairs,1),:) = pairs;
    G.nnc.T(end+1:end+size(pairs,1)) = prod(tpairs,2)/sum([sum(T1),sum(T2)]);
    temp_nnctype = repmat({'star-delta'},size(pairs,1),1);
    G.nnc.type = cat(1,G.nnc.type,temp_nnctype);
end
end
return
