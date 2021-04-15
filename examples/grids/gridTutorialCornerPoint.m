%% How to Read, Display, and Manipulate Corner-Point Grids
% In this example we will show several examples how to manipulate and plot
% corner-point data. As an example, we will use a family of simple models
% with a single fault, which can be realized in different resolutions.


%% Creating and visualizing a small realization of the model
% We start by generating an input stream in the Eclipse format (GRDECL) for
% a 2x1x2 realization of the model
grdecl = simpleGrdecl([2, 1, 2], 0.15)  %#ok (intentional display)

%%
% From the output of  <matlab:help('simpleGrdecl') simpleGrdecl>, we see
% that the file contains four fields:
%
% * The dimension of the underlying logical Cartesian grid (keyword
% SPECGRID, equal 2x1x2)
% * The coordinates of the pillars (keyword COORD, 6x3x2 values)
% * The coordinates along the pillars (keyword ZCORN, 8x2x1x2 values)
% * The flag for active/inactive cells (keyword ACTNUM, 2x1x2 values)
%
% The next step is to process this input stream to build the structure for
% the grid, which here will be represented in a fully unstructured format.
%
% To simplify the processing, a single horizontal layer of artifical cells
% is added one length unit above the top and below the bottom of the model,
% but not touching the model. (This is what is reported in the first line
% of the output). In the following, we therefore work with a 2x1x4 model,
% as shown in the folowing figure
%
% <<topBottomAdded.png>>
%
% Process the grid
G = processGRDECL(grdecl, 'Verbose', true); clear grdecl;

%%
% In the first phase, we process all faces with normals in the logical
% i-direction. Altogether, there should be 3x1x4=12 faces. As we see from
% the figure, there is fault along the second pair of pillars in the
% i-direction. Therefore, only 8 of the 12 faces are reported as regular
% faces. In the next phase of the processing, we treat all pillars pairs
% (fault stacks) that contain at least one set of non-matching
% z-coordinates. Here, there is one such stack and two cells that have a
% non-matching connection. After having split the non-matching faces, the
% stack contains 7 faces (including one face for each of the artificial
% layers) that are added to the list of faces. Then the process is repeated
% in the logical j-direction, in which there are only matching faces.
%
% The processing assumes that there are no faults in the logical
% k-direction and therefore processes only regular connections. In absence
% of inactive or pinched cells, there should be (2+5)x1x2=14 faces.
% Finally, the four artifical cells are removed and the routine correctly
% reports that there are no inactive cells.
%
% The result of the grid processing is a new structure G, outlined
% below
G     %#ok  (intentional display)

%%
% More information about the grid structure can be found in the
% <matlab:help('grid_structure') documentation>.
%
% After the successful processing, we plot the outline of the model and
% mark the faulted faces
clf, subplot('position',[0.025 0.025 0.95 0.95])
plotGrid(G,'FaceColor','b','FaceAlpha', 0.1);
plotFaces(G,find(G.faces.tag>0),'FaceColor','red');
axis equal off, view(40,15)

%% Setting inactive cells
% In this section, we will work with a slightly larger model realization
% and show how one can manipulate the model by setting inactive cells. To
% this end, we start by creating the input stream and building the grid
% structure
grdecl = simpleGrdecl([30, 20, 5], 0.12);
G = processGRDECL(grdecl);

%%
% The resulting model has a single fault and layers given by sinusoidal
% surfaces. First, we plot the outline of the model and mark all fault
% faces in red color
cla
plotGrid(G,'FaceColor','none','EdgeColor',[0.65 0.65 0.65]);
h = plotFaces(G,find(G.faces.tag>0),'FaceColor','red');
axis tight off, view(65,40)

%%
% Then we restrict the model to be inside an ellipse in the logical (i,j)
% coordinates with the principal axis in the i-direction. The model is
% restricted by setting all cells with centerpoint in (i,j)-coordinates
% outside the ellipsis as inactive.
[j,i]  = meshgrid(linspace(-0.95,0.95,G.cartDims(2)), ...
                  linspace(-0.95,0.95,G.cartDims(1)) );
actnum = ones(G.cartDims(1:2));
actnum( (i.^2/4 + j.^2)>0.5) = 0;
actnum = reshape(repmat(actnum,[1 1 G.cartDims(3)]),[],1); clear i j;

set(h, 'FaceAlpha', 0.3);
plotGrid(G,find( actnum(G.cells.indexMap)),'FaceColor','y');
view(20,30), clear h

%%
% To get rid of the inactive cells, we process the model again and plot the
% result to see that it is correct
grdecl.ACTNUM = actnum; clear actnum
G = processGRDECL(grdecl); clear grdecl
cla
plotGrid(G);
axis tight off, view(20,30)

%% Visualizing the layered structure
% To better visualize the layered structure of the model, we plot the
% scalar field 'val(i,j,k)=k' over the grid. This field is constructed
% using the logical Cartesian dimensions of the model, and then the values
% corresponding to the active cells are extracted using 'G.cells.indexMap'.
cla,
val = ones(G.cartDims(1:3));
val = cumsum(val,3);
plotCellData(G,val(G.cells.indexMap),'EdgeColor','k');
clear val

%% Visualizing individual cells
% As our first example, we extract and visualize all cells that are
% neighbor to a fault. To this end, we pick all rows in 'G.faces.neighbors'
% that correspond to faces that are marked as faults, i.e., has a positive
% value for 'G.faces.tag'. If a cell has more than one fault faces, this
% cell will appear more then once in the list. Moreover, outer fault faces
% will be characterized by one of the cell numbers being zero. Before we
% can plot the faces, we must therefore remove redundancy and zero cell
% numbers.
cellList = G.faces.neighbors(G.faces.tag>0, :);
cells    = unique(cellList(cellList>0));
cla,
plotGrid(G,'FaceColor','none','EdgeColor',[0.65 0.65 0.65]);
plotGrid(G,cells);


%%
% Similarly, we can also plot all cells that are neighbors to *internal*
% fault faces by first removing all rows in cellList that have one zero
% entry
cellList = cellList(all(cellList>0,2),:);
cells    = unique(cellList(cellList>0));
plotGrid(G,cells,'FaceColor','g');
clear cellList cells

%% Using ijk-subscripts to access cells
% In the second example, we demonstrate how one can use the logical
% structure of the corner-point grid to pick all active cells corresponding
% to the logical index map [1:2:end,:,1:end/2] and color them using the
% i-value. To be able to access the grid using (i,j,k) adressing, we must
% first construct the subscript values 'ijk' from the linear index
% 'G.cells.indexMap' using the builtin function 'ind2sub'.

% Plot grid outline
cla
plotGrid(G,'FaceColor','none','EdgeColor',[0.65 0.65 0.65]);

% Compute the subindex
clear ijk
[ijk{1:3}] = ind2sub(G.cartDims, G.cells.indexMap);
ijk = [ijk{:}];

% Exctract all active cells for subscript values [1:2:end,:,1:end/2]
[J,K]=meshgrid(1:G.cartDims(2),1:G.cartDims(3)/2);
col = ['b','r','g','c','m','y','w'];
for i=1:2:G.cartDims(1)
   cellNo = find(ismember(ijk, [repmat(i,[numel(J),1]), J(:), K(:)],'rows'));
   plotGrid(G,cellNo,'FaceColor',col(mod(i,numel(col))+1));
end
%clear J K cellNo ijk col

%% Creating a coarse partitioning
% Base on the fine-grid model, we will now create a 5x4x3 coarse
% partitioning with an uneven partitioning in the vertical direction. To
% this end, we will use functionality from the |coarsegrid| module.
mrstModule add coarsegrid
nz = G.cartDims(3);
blockIx = partitionLayers(G,[5 4], ceil([1 nz/2 nz nz+1]));
cla
plotCellData(G,mod(blockIx,8),'EdgeColor','k');
outlineCoarseGrid(G,blockIx, 'FaceColor', 'none', 'LineWidth', 2);
axis tight off, view(20,30)

%%
% Because the partitioning has been performed in logical index space, we
% have so far disregarded the fact the some of the blocks may contain
% disconnected cells because of erosion, faults, etc. Here we see that the
% coarse blocks with indices [3,:,:] do not consist of a single set of
% connected cells. We therefore postprocess the grid to split disconnected
% blocks in two.
blockIx = processPartition(G,blockIx);
cla
plotCellData(G,mod(blockIx,8),'EdgeColor','k');
outlineCoarseGrid(G,blockIx, 'FaceColor', 'none', 'LineWidth', 2);
axis tight off, view(20,30)

%%
% From the plot above, it is not easy to see the shape of the individual
% coarse blocks. We therefore first highlight three representative blocks.
blocks = [4, 6, 28];
col = ['b', 'r', 'g', 'c', 'm', 'y'];

cla
plotGrid(G, 'EdgeColor', [0.75, 0.75, 0.75], 'FaceColor', 'w')
outlineCoarseGrid(G, blockIx, 'FaceColor', 'none', 'LineWidth', 2)

for i = 1 : numel(blocks),
   plotGrid(G, blockIx == blocks(i), ...
            'FaceColor', col(mod(i, numel(col)) + 1));
end

%%
% Alternatively, we can view the partitioning using |explosionView|, which
% moves all blocks out a certain distance from the centerpoint of the
% model.
cla
G = computeGeometry(G);
explosionView(G,blockIx,.5), view(30,50); axis tight off
colormap(colorcube(max(blockIx)))

%% Build the coarse-grid
% Having obtained a partition we are satisfied with, we build the
% coarse-grid structure. This structure consists of three parts:
%
% * the cell structure giving the number of blocks and the indices of the
% cells contained in each block
% * the face structure giving the number of coarse faces and the indices of
% the neighbouring blocks
% * a cellFaces array as in the fine-grid structure
CG = generateCoarseGrid(G, blockIx);
CG             %#ok  (intentional display)
CG.cells       %  (intentional display)
CG.faces       %  (intentional display)

%%
% Let us now use CG to inspect some of the blocks in the coarse grid. To
% this end, we arbitrarily pick a few blocks and inspect these block and
% their neighbours. For the first block, we plot the cells but do not
% highlight the faces that have been marked as lying on a fault
cla
plotBlockAndNeighbors(CG, blocks(1), 'PlotFaults', false([2, 1]))
view(30, 50), axis tight, colormap(jet)

%%
% The next block lies at a fault, so therefore we highlight the fault faces
% as well
cla
plotBlockAndNeighbors(CG, blocks(3), 'PlotFaults', true([2, 1]))
view(20, 30)

%%
displayEndOfDemoMessage(mfilename)

% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
