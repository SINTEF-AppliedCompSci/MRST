%% Example of Realistic Corner-Point Geometry: SAIGUP
% The <http://www.fault-analysis-group.ucd.ie/Projects/SAIGUP.html>SAIGUP
% project is a systematic assessment of uncertainty in reserves
% and production estimates within an objectively defined geological
% parameterisation encompassing the majority of European clastic oil
% reservoirs. A broad suite of shallow marine sedimentological reservoir
% types are indexed to continuously varying 3D anisotropy and heterogeneity
% levels. Structural complexity ranges from unfaulted to compartmentalised,
% and fault types from transmissive to sealing. Several geostatistical
% realisations each for the geologically diverse reservoir types covering
% the pre-defined parameter-space are up-scaled, faulted and simulated with
% an appropriate production strategy for an approximately 20 year period.
% Herein, we will inspect in detail one instance of the model, which can be
% downloaded from the <http://www.sintef.no/Projectweb/MRST>MRST webpage

%% Load data and convert units
% We start by reading the model from a file in the Eclipse format (GRDECL)
% that can be read using the <matlab:doc('readGRDECL') readGRDECL>
% function. (If the model is not available the
% <matlab:doc('getDatasetPath') getDatasetPath> function will download and
% install it for you).
grdecl = fullfile(getDatasetPath('SAIGUP'), 'SAIGUP.GRDECL');
grdecl = readGRDECL(grdecl);
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

%% Inspect full model
G = processGRDECL(grdecl, 'Verbose', true);

%%
G = computeGeometry(G);
rock = grdecl2Rock(grdecl, G.cells.indexMap);

%% Inspect geometry
% Construct pillars and corner-points to verify the resolution of the model
[X,Y,Z] = buildCornerPtPillars(grdecl,'Scale',true);            %#ok<NASGU>
dx = unique(diff(X)).'                                          %#ok<NOPTS>
[x,y,z] = buildCornerPtNodes(grdecl);                           %#ok<ASGLU>
dz = unique(reshape(diff(z,1,3),1,[]))                          %#ok<NOPTS>
clear x y z X Y Z

%% 
% Next, we look at small faces created to enforce a matching grid. Vertical
% faces that are not subdividided will be parallelograms with a 300 m^2
% area. We therefore look at faces having smaller area
hist(G.faces.areas(G.faces.areas<300),100);
xlabel('Face area, m^2');
disp(['Area less than 0.01 m^2: ', num2str(sum(G.faces.areas<0.01))])
disp(['Area less than 0.10 m^2: ', num2str(sum(G.faces.areas<0.1))])
disp(['Area less than 1.00 m^2: ', num2str(sum(G.faces.areas<1))])

%%
% Plot where these faces appear in the model
clf
plotGrid(G,'FaceColor','none','EdgeAlpha',0.1);
plotFaces(G,find(G.faces.tag>0 & G.faces.areas>=290),'y','edgea',0.1);
plotFaces(G,find(G.faces.areas<290),'r','edgea',0.1);
view(-105,45), axis tight off

%%
% Process the model using a tolerance to get rid of the smallest faces
for a=[0.05 0.1 .25 0.5 1]
   G2 = processGRDECL(grdecl, 'Tolerance', a);
   G2 = computeGeometry(G2);
   disp([ a, min(G2.faces.areas)])
end
clear G2


%% The layered structure
% Show layered structure using a very simple trick: create a matrix with
% ones in all cells of the logical Cartesian grid and then do a cummulative
% summation in the vertical direction to get increasing values.
clf, args = {'EdgeColor','k','EdgeAlpha',0.1};
val = cumsum(ones(G.cartDims),3);
plotCellData(G,val(G.cells.indexMap),args{:});
view(-80,50), axis off tight

%%
% Unfortunately, our attempt at visualizing the layered structure was not
% very successful. We therefore try to extract and visualize only the cells
% that are adjacent to a fault.
cellList = G.faces.neighbors(G.faces.tag>0, :);
cells    = unique(cellList(cellList>0));
cla,
plotCellData(G,val(G.cells.indexMap),cells,args{:});
view(-120,40)

%%
% Restrict the plot to only parts of the grid using ismember, which uses
% O(n log n) operations
clear ijk
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap); ijk = [ijk{:}];
[I,J,K] = meshgrid(1:9,1:30,1:20);
bndBox = find(ismember(ijk,[I(:), J(:), K(:)],'rows'));
inspect = cells(ismember(cells,bndBox));
cla,
plotCellData(G,val(G.cells.indexMap), inspect,'EdgeColor','k');
axis tight off

%%
% Restrict the plot to only parts of the grid using logical indexing which
% uses O(n) operations
clear ijk
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
I = false(G.cartDims(1),1); I(1:9)=true;
J = false(G.cartDims(2),1); J(1:30)=true;
K = false(G.cartDims(3),1); K(1:20)=true;

pick = I(ijk{1}) & J(ijk{2}) & K(ijk{3});
pick2 = false(G.cells.num,1); pick2(cells) = true;
inspect = find(pick & pick2);

cla,
plotCellData(G,val(G.cells.indexMap), inspect,'EdgeColor','k');
axis tight off

%%
% Let us now try to also colour the faces of the model
cellno  = gridCellNo(G);
faces   = unique(G.cells.faces(pick(cellno), 1));
inspect = faces(G.faces.tag(faces)>0);
plotFaces(G,inspect,[.7 .7 .7],'edgec','r');

%% Making a ribbon plot
% 
clear ijk;
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);

% Ribbon in x-direction
I = false(G.cartDims(1),1); I(1:5:end)=true;
J = true(G.cartDims(2),1);
K = true(G.cartDims(3),1);

pickX = I(ijk{1}) & J(ijk{2}) & K(ijk{3});

% Ribbon in y-direction
I = true(G.cartDims(1),1);
J = false(G.cartDims(2),1); J(1:10:end) = true;
K = true(G.cartDims(3),1);

pickY = I(ijk{1}) & J(ijk{2}) & K(ijk{3});

clf,
plotCellData(G,rock.poro, pickX | pickY,'EdgeColor','k','EdgeAlpha',0.1);
view(-100,45), axis tight off

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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
