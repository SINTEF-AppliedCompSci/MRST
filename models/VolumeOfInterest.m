classdef VolumeOfInterest
    % Class for Volume of Interest (VOI) in the Corner-point grid or Cartesian grid
    properties
        CPG            % Corner-point or Cartesian grid struct
        well           % Struct for well information
        boundary       % 2D VOI boundary specified by the polygon
        extraLayers    % Extra layers above and below the layers where the
                       % well is located
    end
    
    methods
        
        function volume = VolumeOfInterest(G, well, pbdy, nextra, varargin)
            volume.CPG         = G;
            volume.well        = well;
            volume.boundary    = pbdy;
            % nextra(1): Layer numbers above the layers that well occupies
            % nextra(2): Layer numbers below the layers that well occupies
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
            % Get logical indices of Grid. The Grid should be structured in
            % all three dimensions.
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
            % Find face indicator at layered dimension (typically is 'Z')
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
        
        function wc = PeacemanWellCells(volume, varargin)
            % Find well cells of Peaceman's well model, require 'wellpaths'
            % module
            pW   = volume.well.trajectory;
            wph  = makeSingleWellpath(pW);
            wc   = findWellPathCells(volume.CPG, wph);
        end
        
        function ij = ijIndicesFromBoundary(volume, varargin)
            % Get i and j indices from the defined 2D boundary
            pbdy = volume.boundary;
            % Get ij indices of the volume from a defined 2D polygon
            % All well points should be located inside the boundary
            pW = volume.well.trajectory;
            assert(all(inpolygon(pW(:,1), pW(:,2), pbdy(:,1), pbdy(:,2))), ...
                ['Well points outside the boundary were defected, ', ...
                'try to enlarge the boundary']);
            [I, J, K] = volume.logicalIndices();
            % Find VOI cells per layers
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
            tabi = tabulate(ij(:,1));
            badi = ismember(ij(:,1), tabi(tabi(:,2) == 1, 1));
            tabj = tabulate(ij(:,2));
            badj = ismember(ij(:,2), tabj(tabj(:,2) == 1, 1));
            ij = ij(~(badi|badj), :);
        end
        
        function k = kIndicesFromExtraLayers(volume, varargin)
            % Get the grid layer index from the extra layers and layers
            % occupied by the well
            nex = volume.extraLayers;
            % Get k from grid layers that well occupies and extra defined
            % layers. nextra = number of layers [above, below] the well
            % layers
            [~, ~, K] = volume.logicalIndices();
            wc    = volume.PeacemanWellCells();
            kwc   = K(wc);
            kmin  = min(kwc) - nex(1);
            kmax  = max(kwc) + nex(2);
            k = (kmin : kmax)';
            k = k( k <= max(K) & k >= min(K) );
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
            % For Corner-point grid, the layer-faces are Z- and Z+ faces
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
            % Get all boundary nodes and faces of the volume
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
            % Boundary faces of g, sorted, counter clock wise
            bfg = find( ~all(Ng,2) );
            bfg = sortPtsCounterClockWise(g.faces.centroids, {bfg});
            bfg = bfg{1};
            % Nodes of bf, also boundary nodes of g
            [bfng, pos] = gridFaceNodes(g, bfg);
            assert(all(diff(pos)==2))
            bfng = reshape(bfng, 2, [])';
            % Boundary nodes of VOI in g, sorted, counter clock wise
            bng = arrayfun(@(r)bfng(r, ~ismember(bfng(r,:), bfng(r-1,:))), ...
                (2:size(bfng,1)-1)');
            idx  = ismember(bfng(1,:), bfng(2,:));
            bng  = [bfng(1,~idx); bfng(1,idx); bng];
            % Boundary nodes of VOI in G, sorted, counter clock wise
            bn  = nu(bng);
            % The following are the preparations of building pebi grid
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
            % Insert the Z cells
            bcg = [bcg, cc]';
            bcg = bcg(:);
            bcg = bcg(bcg~=0);
            % Boundary faces of VOI in g,  cell index of g equals to face
            % index of VOI f
            bf = f(bcg);
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
            % Plot layer faces of the volume
            figure, hold on, axis off
            arrayfun(@(x)plotFaces(volume.CPG, packed.faces{x}, 'facecolor',...
                rand(3,1)), (1:numel(packed.faces))')
            title('Layer faces of the volume')
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
            %  WR: The 2D well region that consists of should include the
            %  following field:
            %   'ny' - The number of Cartesian cells in Y direction
            %   'ny' - The size of Cartesian region in Y direction
            %   'na' - The number of angular cells in radial region
            %
            % RETURNS:
            % The expanded WR with node and topology information
            %  'points'    - 2D coordinates WR nodes
            %  'connlist'  - Connectivity list of the WR nodes
            %  'connlistC' - Connectivity list of the WR nodes,
            %                only include the Cartesian region
            %  'bdnodes'   - Indices for boundary nodes of WR
            %  'bdnodesC'  - Indices for boundary nodes of WR,
            %                only include the Cartesian region
            %  'cartDims'  - Dimensions of the Cartesian region, [nx, ny]
            %
            pW = volume.well.trajectory;
            nx = volume.well.segmentNum;
            ny = WR.ny;
            if mod(ny,2) == 1
                warning(['ny must be an even number, ny+1 (%d) is',...
                    'used instead'], ny+1)
                ny = ny + 1;
            end
            ly = WR.ly;
            na = WR.na;
            % Get points corresponding to single well node
            p = arrayfun(@(ii)pointsSingleWellNode(pW, ly, ny, na, ii), ...
                (1:nx+1)');
            
            pall = [vertcat(p.cart); vertcat(p.rad)];
            pbdy = volume.boundary;
            assert(all(inpolygon(pall(:,1), pall(:,2), pbdy(:,1), pbdy(:,2))),...
                ['Points outside the boundary were detected, please reduce ',...
                'the size of Cartesian region']);
            [t, tC, bn, bnC] = getConnListAndBdyNodeWR2D(p, ny, na);
            
            % Asssign data to WR
            WR.points     = p;
            WR.connlist   = t;
            WR.connlistC  = tC;
            WR.bdnodes    = bn;
            WR.bdnodesC   = bnC;
            WR.cartDims   = [nx, ny];
        end
        
        function GV = ReConstructToUnstructuredGrid(volume, WR, layerRf, varargin)
            % Reconstruct CPG in VOI to layered unstructured grid.
            % We use the open-source triangle generator 'DistMesh'
            % (Per-Olof Persson) to obtain high-quality triangles.
            % We adopt the scaled edge length function:
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
            %  layerRf - Number of refined layers in each VOI layer
            %
            %  Optional:
            %   'multiplier'- Multiplier of the scaled edge length function
            %   'maxIter'   - The maximum number of distmesh iterations
            %   'gridType'  - Grid type, 'Voronoi' (default) | 'triangular'
            %
            % RETURNS:
            %  GV - Layered unstructured VOI grid
            %
            fprintf(' -- Reconstructing the CPG to unstructured VOI grid\n')
            
            opt = struct('multiplier', 0.2,...
                'maxIter', 500,...
                'gridType', 'triangular');
            opt = merge_options(opt, varargin{:});
            
            packed = volume.allInfoOfVolume();
            if ~isfield(WR, 'points')
                WR = volume.prepareWellRegionNodes2D(WR);
            end
            
            % Generate nodes for the unstructured grid
            [players, t, bdyID] = ...
                generateVOIGridNodes(volume.CPG, packed, WR, layerRf, opt);
            
            % Construct the 2D grid
            p = players{1}(:, [1,2]);
            t = sortPtsCounterClockWise(p, t);
            gV = tessellationGrid(p, t);
            gV.nodes.boundary = bdyID;
            gV.cartDims = WR.cartDims;
            
            % Extrude the 2D grid to 3D grid
            GV = makeLayeredGridNWM(gV, players, 'connectivity', t);
            if size(layerRf,1) < size(layerRf,2)
                layerRf = layerRf';
            end
            GV.layers.refinement = layerRf(1:numel(packed.cells));
            GV.parentInfo = packed;
            
            dispInfo(GV);
        end
        
        function plot2DWRSubGrid(volume, WR)
            % Plot the subgrid of2D well region
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
            % Compute maximum lengths of well segments
            pW = volume.well.trajectory(:, [1,2]);
            L = diff(pW, 1);
            L = sqrt(sum(L.^2,2));
            fprintf('    Info : The maximum well-segment length in 2D is %.2f\n', max(L))
        end
        
    end
end

function  plot3DLines(p, mk, lw)
    plot3(p(:,1), p(:,2), p(:,3), mk, 'linewidth', lw);
end


