function Gl = frac_frac_nnc_extruded(G, Gl, F, fracture, flayers)
% frac_frac_nnc_extruded assigns NNC connections to fracture-fracture
% intersections (similar to its 2D counterpart frac_frac_nnc) and also
% assigns a transmissibility to each connection using the star-delta
% transformation. See see SPE-88812-PA, Karimi-Fard et al, 2004. The NNC's
% are first identified in 2D and then extruded to 3D.

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
    findex2D = []; findex3D = [];
    for i = 1:numel(F)
        findex2D = [findex2D,F(i).cells.start]; %#ok
        findex3D = [findex3D, ...
            Gl.FracGrid.(['Frac',num2str(i)]).cells.start]; %#ok
    end
    findex2D(end+1) = findex2D(end) + F(i).cells.num;
    findex3D(end+1) = findex3D(end) + Gl.FracGrid.(['Frac',num2str(i)]).cells.num;
    for i = 1:size(fracture.intersections.lines,1) %% 2D for now
        lines = fracture.intersections.lines(i,:);
        coords = fracture.intersections.coords(i,:);
        Gf1 = G.FracGrid.(['Frac',num2str(lines(1))]);
        Gf2 = G.FracGrid.(['Frac',num2str(lines(2))]);
        %
        [~,Gface(1)] = ismember(roundsd(coords,5),roundsd(F(lines(1)).nodes.coords,5),'rows');
        [~,Gface(2)] = ismember(roundsd(coords,5),roundsd(F(lines(2)).nodes.coords,5),'rows');
        %  
        cells_l1 = cell(2,1);
        for j = 1:numel(Gface)
            if Gface(j) == 1
                cells_l1{j,1} = Gface(j);
            elseif Gface(j) == size(F(lines(j)).nodes.coords,1)
                cells_l1{j,1} = Gface(j)-1;
            else
                cells_l1{j,1} = [Gface(j)-1,Gface(j)];
            end
        end
        cells_l1{1,1} = Gf1.cells.start-1 + cells_l1{1,1};
        cells_l1{2,1} = Gf2.cells.start-1 + cells_l1{2,1};
        diff1 = diff([F(lines(1)).nodes.coords(1,:);F(lines(1)).nodes.coords(end,:)]);
        diff2 = diff([F(lines(2)).nodes.coords(1,:);F(lines(2)).nodes.coords(end,:)]);
        if diff1(1) == 0
            Gface(1) = Gf1.faces.num - (size(F(lines(1)).nodes.coords,1)-Gface(1));
        end
        if diff2(1) == 0
            Gface(2) = Gf2.faces.num - (size(F(lines(2)).nodes.coords,1)-Gface(2));
        end
        %
        cent1 = Gf1.faces.centroids(Gface(1),:);
        cent2 = Gf2.faces.centroids(Gface(2),:);
        %
        Gf1 = Gl.FracGrid.(['Frac',num2str(lines(1))]);
        Gf2 = Gl.FracGrid.(['Frac',num2str(lines(2))]);
        %
        faces1 = find(ismember( ...
           roundsd(Gf1.faces.centroids(:,1:2),12),...
           roundsd(cent1,12),'rows'));
        faces2 = find(ismember(...
           roundsd(Gf2.faces.centroids(:,1:2),12),...
           roundsd(cent2,12),'rows'));
        %-----------Get half face transmissibilities----------------------%
        Trans1 = computeTrans(Gf1,Gf1.rock);
        Trans2 = computeTrans(Gf2,Gf2.rock);
        
        for j = 1:numel(flayers)
            cells_l = cell(2,1);
            
            ind2D = cells_l1{1,1} - findex2D(lines(1));
            ind3D = findex3D(lines(1)) + ind2D;
            cells_l{1,1} = ind3D + (j-1)*F(lines(1)).cells.num;
            
            ind2D = cells_l1{2,1} - findex2D(lines(2));
            ind3D = findex3D(lines(2)) + ind2D;
            cells_l{2,1} = ind3D + (j-1)*F(lines(2)).cells.num;
            

            T1 = Trans1(Gf1.cells.faces(:,1) == faces1(j)); % Half-trans multiplers. t1(1) => for cells_l1(1)
            T2 = Trans2(Gf2.cells.faces(:,1) == faces2(j)); % Half-trans multiplers. t2(1) => for cells_l2(1)
            
            %----Define NNC's and respective trans using star-delta-------%
            
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
            Gl.nnc.cells(end+1:end+size(pairs,1),:) = pairs;
            Gl.nnc.T(end+1:end+size(pairs,1)) = prod(tpairs,2)/sum([sum(T1),sum(T2)]);
            temp_nnctype = repmat({'star-delta'},size(pairs,1),1);
            Gl.nnc.type = cat(1,Gl.nnc.type,temp_nnctype);
        end
    end
end
return
