classdef VolumeOfInterest
% Class for volume of interest (VOI) in the Corner-point grid (CPG) or 
% Cartesian grid which generates the geometrical information of CPG or 
% Cartesian grid in VOI and constructs the unstructured VOI grid.
    
    properties
        CPG            % CPG or Cartesian grid structure
        well           % Structure of well information
        boundary       % 2D VOI boundary specified by the polygon
        extraLayers    % Extra layers above and below the layers where the
                       % well is located
    end
    
    methods
        
        function volume = VolumeOfInterest(G, well, pbdy, nextra, varargin)
            volume.CPG         = G;
            volume.well        = well;
            volume.boundary    = pbdy;
            % nextra(1): Number of layers above the HW layers
            % nextra(2): Number of layers below the HW layers
            volume.extraLayers = nextra;
            volume.plotVolumeBoundaries(1, 'plotClippedBoundary',false)
            % All well points should be located inside the boundary
            pW = volume.well.trajectory;
            if ~all(inpolygon(pW(:,1), pW(:,2), pbdy(:,1), pbdy(:,2)))
                error(['Well points outside the boundary are ', ...
                    'defected, please enlarge the boundary']);
            end
        end
        
        function varargout = logicalIndices(volume, varargin)
            % Get logical indices of CPG
            G = volume.CPG;
            if nargin > 1
                indices = gridLogicalIndices(G, varargin{1});
            else
                indices = gridLogicalIndices(G);
            end
            if nargout == 1
                varargout = indices;
            else
                varargout = arrayfunUniOut(@(x)indices{x}, (1:nargout)');
            end
        end
        
        function indicator = layerFaceIndicator(volume, varargin)
            % Find face indicator of layered dimension (typically is 'Z')
            G = volume.CPG;
            if isfield(G, 'cartDims') && size(G.cells.faces,2)==2
                indicator = [5, 6];
            else
                indicator = [nan, nan];
            end
        end
        
        function c = logicalToArray(volume, ijk, varargin)
            % Convert logical indices to array indices
            [I, J, K] = volume.logicalIndices();
            c = arrayfun(@(x)find(I == ijk(x,1) & J == ijk(x,2) & ...
                K == ijk(x,3)), (1:size(ijk,1))');
        end
        
        function packed = allInfoOfVolume(volume, varargin)
            % Get all information of the volume
            packed.cells = volume.cellsOfVolume();
            packed.faces = volume.facesOfVolume(packed.cells);
            packed.nodes = volume.nodesOfVolume(packed.faces);
            [bn, bf] = volume.boundaryInfoOfVolume(packed.faces);
            packed.bdyNodes = bn;
            packed.bdyFaces = bf;
            packed.boxCells = volume.boxCellsOfVolume();
            packed.boxFaces = volume.facesOfVolume(packed.boxCells);
            packed.PeacemanCells   = volume.PeacemanWellCells();
            packed.clippedBoundary = cellfunUniOut(...
                @(n)volume.CPG.nodes.coords(n, [1,2]), bn);
            packed.KIndices = volume.kIndicesFromExtraLayers();
        end
        
        function wc = PeacemanWellCells(volume, varargin)
            % Find well cells of Peaceman's well model, require 'wellpaths'
            % module
            pW   = volume.well.trajectory;
            wph  = makeSingleWellpath(pW);
            wc   = findWellPathCells(volume.CPG, wph);
        end
        
        function ij = ijIndicesFromBoundary(volume, varargin)
            % Get i and j indices f the volume from the defined 2D boundary
            pbdy = volume.boundary;
            % All well points should be located inside the boundary
            pW = volume.well.trajectory;
            assert(all(inpolygon(pW(:,1), pW(:,2), pbdy(:,1), pbdy(:,2))), ...
                ['Well points outside the boundary were defected, ', ...
                'try to enlarge the boundary']);
            [I, J, K] = volume.logicalIndices();
            % Find VOI cells per layer
            c = cell(max(K), 1);
            for k = min(K) : max(K)
                cktol = find(K == k);
                xy    = volume.CPG.cells.centroids(cktol, [1, 2]);
                in    = inpolygon(xy(:,1), xy(:,2), pbdy(:,1), pbdy(:,2));
                c{k}  = cktol(in);
            end
            % Combine the cells and extract the logical indices
            ij = cellfunUniOut(@(c)[I(c), J(c)], c);
            ij = cell2mat(ij);
            ij = unique(ij, 'rows');
            % Remove 'bad' ij (appears only once)
            tabi = tabulate_NWM(ij(:,1));
            badi = ismember(ij(:,1), tabi(tabi(:,2) == 1, 1));
            tabj = tabulate_NWM(ij(:,2));
            badj = ismember(ij(:,2), tabj(tabj(:,2) == 1, 1));
            ij = ij(~(badi|badj), :);
        end
        
        function k = kIndicesFromExtraLayers(volume, varargin)
            % Get the grid layer indices from the extra layers and layers
            % occupied by the well
            % nex = number of layers [above, below] the well layers
            nex = volume.extraLayers;
            [~, ~, K] = volume.logicalIndices();
            wc    = volume.PeacemanWellCells();
            kwc   = K(wc);
            kmin  = min(kwc) - nex(1);
            kmax  = max(kwc) + nex(2);
            k = (kmin : kmax)';
            k = k( k <= max(K) & k >= min(K) );
        end
        
        function c = cellsOfVolume(volume, varargin)
            % Get all cells of the volume
            k = volume.kIndicesFromExtraLayers();
            c = arrayfunUniOut(@(k)volume.getCellsSingleLayer(k), k);
        end
        
        function c = getCellsSingleLayer(volume, k, varargin)
            % Get layer-k cells inside the defined 2D polygon
            ij  = volume.ijIndicesFromBoundary();
            ijk = [ij, k*ones(size(ij,1),1)];
            c   = volume.logicalToArray(ijk);
        end
        
        function f = facesOfVolume(volume, c, varargin)
            % Get all layer-faces of the volume
            indicator = volume.layerFaceIndicator();
            f = cellfunUniOut(@(c)volume.getLayerFacesFromCells...
                (c, indicator(1)), c);
            f0 = volume.getLayerFacesFromCells(c{end},indicator(2));
            f  = [f; f0];
        end
        
        function f = getLayerFacesFromCells(volume, c, indicator, varargin)
            % Get layer-faces of cell c in single layer
            % For CPG or Cartesian grid, the layer-faces are Z- and Z+ 
            % faces
            G   = volume.CPG;
            cf  = G.cells.faces(mcolon(G.cells.facePos(c), ...
                G.cells.facePos(c+1) - 1),1);
            dir = G.cells.faces(mcolon(G.cells.facePos(c), ...
                G.cells.facePos(c+1) - 1),2);
            f   = cf(dir == indicator);
        end
        
        function n = nodesOfVolume(volume, f, varargin)
            % Get all nodes of layer-faces of the volume
            n = cellfunUniOut(@(f)volume.getNodesFromFaces(f), f);
        end
        
        function n = getNodesFromFaces(volume, f)
            % Get nodes of layer-faces on single surface
            n = arrayfunUniOut(@(f)gridFaceNodes(volume.CPG, f), f);
        end
        
        function boxc = boxCellsOfVolume(volume, varargin)
            % Get all box cells of the volume
            k = volume.kIndicesFromExtraLayers();
            boxc = arrayfunUniOut(@(k)volume.getBoxCellsSingleLayer(k), k);
        end
        
        function boxc = getBoxCellsSingleLayer(volume, k, varargin)
            % Get layer-k box cells. The defined 2D boundary located inside
            % the box. This will be useful in the interpolations of rock
            % properties.
            [I, J, K] = volume.logicalIndices();
            ij  = volume.ijIndicesFromBoundary();
            en  = 3;
            [imin, imax, jmin, jmax] = deal(min(ij(:,1)) - en, ...
                max(ij(:,1)) + en,...
                min(ij(:,2)) - en,...
                max(ij(:,2)) + en);
            boxc = find(I>=imin & I<=imax & J>=jmin & J<=jmax & K==k);
        end
        
        function [bn, bf] = boundaryInfoOfVolume(volume, f, varargin)
            % Get all boundary nodes and layer-faces of the volume
            [bn, bf] = cellfun(@(f)volume.getBoundaryInfoSingleSurface(f),...
                f, 'UniformOutput', false);
        end
        
        function [bn, bf] = getBoundaryInfoSingleSurface(volume, f, varargin)
            % Get boundary information of faces on single surface
            % bn:  sorted boundary nodes
            % bf:  sorted boundary faces
            G = volume.CPG;
            n = arrayfunUniOut(@(f)gridFaceNodes(G, f), f);
            n = sortPtsCounterClockWise(G.nodes.coords(:,1:2), n);
            assert(all(cellfun(@numel, n)==4))
            % Build a local grid g to find boundary nodes of G
            nd = cell2mat(n);
            [nu, ~, ic] = unique(nd);
            t = reshape(ic, 4, [])';
            p = G.nodes.coords(nu, [1, 2]);
            g = tessellationGrid(p, t);
            g = computeGeometry(g);
            Ng = g.faces.neighbors;
            % Boundary faces of g, sorted, counterclockwise
            bfg = find( ~all(Ng,2) );
            bfg = sortPtsCounterClockWise(g.faces.centroids, {bfg});
            bfg = bfg{1};
            % Nodes of bf, also boundary nodes of g
            [bfng, pos] = gridFaceNodes(g, bfg);
            assert(all(diff(pos)==2))
            bfng = reshape(bfng, 2, [])';
            % Boundary nodes of VOI in g, sorted, counterclockwise
            bng = arrayfun(@(r)bfng(r, ~ismember(bfng(r,:), bfng(r-1,:))), ...
                (2:size(bfng,1)-1)');
            idx  = ismember(bfng(1,:), bfng(2,:));
            bng  = [bfng(1,~idx); bfng(1,idx); bng];
            % Boundary nodes of VOI in G, sorted, counterclockwise
            bn  = nu(bng);
            % The following are the preparations of building Voronoi grid
            % Boundary cells in g
            bcg = sum(Ng(bfg, :),2);
            bcg = unique(bcg, 'stable');
            N = [bcg, [bcg(2:end); bcg(1)]];
            cc = zeros(size(bcg));
            for ii = 1 : size(N,1)
                c1 = Ng(any(Ng == N(ii,1),2), :);
                c1 = unique(c1); c1 = c1(c1~=0 & c1~=N(ii,1));
                c2 = Ng(any(Ng == N(ii,2),2), :);
                c2 = unique(c2); c2 = c2(c2~=0 & c2~=N(ii,2));
                if ~isempty(intersect(c1, c2))
                    cc(ii) = intersect(c1, c2);
                end
            end
            % Insert the 'Z' cells
            bcg = [bcg, cc]';
            bcg = bcg(:);
            bcg = bcg(bcg~=0);
            % Boundary faces of VOI in g, cell index of g equals to face
            % index of VOI face
            bf = f(bcg);
            assert(numel(bf) == unique(numel(bf)), ['Isolate boundary',...
                ' faces are detected, please redefine the boundary polygon'])
        end
        
        function WR = prepareWellRegionNodes2D(volume, WR)
            % Prepare the 2D well region nodes.
            % The unstructured VOI grid includes a 2D well region (WR). The
            % WR is composed of a Cartesian region and two half-radial
            % regions in xy plane, which are used to connect the HW grid.
            % This function generates the grid nodes for the three
            % structured regions.
            % For the Cartesian region, the X axis extends along the well
            % trajectory:
            %     ---------> X
            %    |    -----------------------------------
            %  Y |    --------- Well trajectory ---------
            %    V    -----------------------------------
            %
            % PARAMETERS:
            % WR - The 2D WR structure that consists of following fields:
            %  'ny' - The number of Cartesian cells in Y direction
            %  'ny' - The size of Cartesian region in Y direction
            %  'na' - The number of angular cells in radial region
            %
            % RETURNS:
            % The expanded WR structure with node and topology information:
            %  'points'    - 2D coordinates WR nodes
            %  'connlist'  - Connectivity list of the whole well region
            %  'connlistC' - Connectivity list of the Cartesian region
            %  'bdnodes'   - Indices of outer boundary nodes of the whole 
            %                well region
            %  'bdnodesC'  - Indices of outer boundary nodes of the 
            %                Cartesian region
            %  'cartDims'  - Dimensions of the Cartesian region, [nx, ny]
           
            pW = volume.well.trajectory;
            nx = volume.well.segmentNum;
            ny = WR.ny;
            for i = 1 : numel(ny)
                if mod(ny(i),2) == 1
                    warning(['ny must be an even number, ny(%d)+1 [%d] is',...
                        ' used instead'], i, ny(i)+1)
                    ny(i) = ny(i) + 1;
                end
            end
            ly = WR.ly;
            na = WR.na;
            
            % Generate the WR points corresponding to all well nodes
            p = arrayfun(@(ii)pointsSingleWellNode(pW, ly, ny, na, ii), ...
                (1:nx+1)');
            pall = [vertcat(p.cart); vertcat(p.rad)];
            pbdy = volume.boundary;
            assert(all(inpolygon(pall(:,1), pall(:,2), pbdy(:,1), pbdy(:,2))),...
                ['Points outside the boundary were detected, please reduce ',...
                'the size of Cartesian region']);
            
            % Get connectivity list and boundary nodes of WR
            [t, tC, bn, bnC] = getConnListAndBdyNodeWR2D(p, sum(ny), na);
            
            % Asssign data to WR
            WR.points     = p;
            WR.connlist   = t;
            WR.connlistC  = tC;
            WR.bdnodes    = bn;
            WR.bdnodesC   = bnC;
            WR.cartDims   = [nx, sum(ny)];
        end
        
        function GV = ReConstructToUnstructuredGrid(volume, WR, layerRf, varargin)
            % Reconstruct CPG in VOI to layered unstructured grid.
            % The open-source triangle generator 'DistMesh' 
            % (Per-Olof Persson) is used to obtain high-quality triangles. 
            % The scaled edge length function is defined as:
            % h(p) = max(multiplier*d(p) +lIB, lOB)
            % to let the point density increases from inner boundary to
            % outer boundary
            % lIB: average length of inner boundary (outer-boundary of WR subgrid)
            % lOB: average length of outer boundary (VOI clipped boundary)
            %
            % PARAMETERS:
            %  WR      - The 2D well region struct, including the node and
            %            topology information. Generated by
            %            'prepareWellRegionNodes2D'
            %  layerRf - Number of refinement layers in each VOI layer
            %
            % KEYWORD ARGUMENTS:
            %   'multiplier'- Multiplier of the scaled edge length function
            %   'maxIter'   - The maximum number of DistMesh iterations
            %   'gridType'  - Grid type, 'Voronoi' (default) | 'triangular'
            %
            % RETURNS:
            %  GV - Layered unstructured VOI grid
           
            fprintf(' -- Reconstructing the CPG to unstructured VOI grid\n')
            
            opt = struct('multiplier', 0.2,...
                'maxIter', 500,...
                'gridType', 'triangular');
            opt = merge_options(opt, varargin{:});
            
            % Geometrical info of VOI
            packed = volume.allInfoOfVolume();
            if ~isfield(WR, 'points')
                WR = volume.prepareWellRegionNodes2D(WR);
            end
            
            % Generate nodes for the unstructured grid
            [pSurfs, t, bdyID] = ...
                generateVOIGridNodes(volume.CPG, packed, WR, layerRf, opt);
            
            % Construct the 2D grid
            p = pSurfs{1}(:, [1,2]);
            t = sortPtsCounterClockWise(p, t);
            gV = tessellationGrid(p, t);
            gV.nodes.boundary = bdyID;
            gV.cartDims = WR.cartDims;
            
            % Extrude the 2D grid to 3D grid
            GV = makeLayeredGridNWM(gV, pSurfs, 'connectivity', t);
            if size(layerRf,1) < size(layerRf,2)
                layerRf = layerRf';
            end
            GV.layers.refinement = layerRf(1:numel(packed.cells));
            GV.parentInfo = packed;
            
            dispInfo(GV);
        end
        
        function plotVolumeBoundaries(volume, packed, varargin)
            % Plot the user defined boundary and clipped boundary
            opt = struct('plotClippedBoundary', true);
            opt = merge_options(opt, varargin{:});
            pW  = volume.well.trajectory;
            pB  = volume.boundary;
            G   = volume.CPG;
            figure, hold on, axis off
            plot3DLines(pW, 's-', 1)
            pBClose = [[pB; pB(1,:)], pW(1,3)*ones(size(pB,1)+1,1)];
            plot3DLines(pBClose, 'o-', 1.5)
            if opt.plotClippedBoundary
                pBClipped = G.nodes.coords(packed.bdyNodes{1},:);
                plot3DLines(pBClipped, 'g.-', 1.5)
                legend({'Well trajectory', 'Specified VOI boundary', ...
                    'Clipped VOI boundary'})
            else
                legend({'Well trajectory', 'Specified VOI boundary'})
            end
            plotGrid(G, 'facecolor', 'none')
        end
        
        function plotVolumeCells(volume, packed)
            % Plot cells inside the volume
            figure, hold on, axis off
            arrayfun(@(x)plotGrid(volume.CPG, packed.cells{x}, 'facecolor',...
                rand(3,1)), (1:numel(packed.cells))')
            title('Cells inside the volume')
        end
        
        function plotVolumeLayerFaces(volume, packed)
            % Plot layer-faces of the volume
            figure, hold on, axis off
            arrayfun(@(x)plotFaces(volume.CPG, packed.faces{x}, 'facecolor',...
                rand(3,1)), (1:numel(packed.faces))')
            title('Layer-faces of the volume')
        end
        
        function plot2DWRSubGrid(volume, WR)
            % Plot the subgrid of the 2D well region
            if ~isfield(WR, 'points')
                WR = volume.prepareWellRegionNodes2D(WR);
            end
            figure, hold on, axis off
            pWR = [vertcat(WR.points.cart); vertcat(WR.points.rad)];
            gWR = tessellationGrid(pWR, WR.connlist);
            plot3DLines(volume.well.trajectory, 'rs-', 2)
            legend('Well trajectory')
            plotGrid(gWR, 'facecolor', 'none')
            title('2D well reigon subgrid')
        end
        
        function maxWellSegLength2D(volume)
            % Display the maximum 2D length of well segments
            pW = volume.well.trajectory(:, [1,2]);
            L = diff(pW, 1);
            L = sqrt(sum(L.^2,2));
            fprintf('    Info : The maximum well-segment length in 2D is %.2f\n', max(L))
        end
        
        function volumeLayerNumber(volume)
            % Display the number of volume layers
            k = volume.kIndicesFromExtraLayers();
            fprintf('    Info : The number of VOI layers is %d (', numel(k))
            arrayfun(@(k)fprintf(' %d ', k), k)
            fprintf(')\n')
        end
        
    end
end

function  plot3DLines(p, mk, lw)
    plot3(p(:,1), p(:,2), p(:,3), mk, 'linewidth', lw);
end


