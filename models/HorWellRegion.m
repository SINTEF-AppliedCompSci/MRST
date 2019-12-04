classdef HorWellRegion
% Class for horizontal well (HW) region in volume of interest (VOI) grid 
% which generates the geometrical information of VOI grid and constructs 
% the radial HW grid
    
    properties
        GVOI             % Layered unstructured VOI grid
        regionIndices    % Logical indices of HW region in VOI grid
        well             % Structure of well information
    end
    
    methods
        
        function hw = HorWellRegion(G, well, regionIndices, varargin)
            % The HW grid is built inside the Cartesian region of VOI grid. 
            % And the logical indices of HW region are specified.
            %
            %      1   ymin                     ymax   ny
            %    ----- ----- ----- ----- ----- ----- -----
            %   |     |     |     |     |     |     |     |    1
            %    ----- ----- ----- ----- ----- ----- -----
            %   |     |  *  |  *  |  *  |  *  |  *  |     |   zmin
            %    ----- ----- ----- ----- ----- ----- -----
            %   |     |  *  |  *  |  *  |  *  |  *  |     |
            %    ----- ----- ----- ----- ----- ----- -----
            %   |     |  *  |  *  |  *  |  *  |  *  |     |   zmax
            %    ----- ----- ----- ----- ----- ----- -----
            %   |     |     |     |     |     |     |     |    nz
            %    ----- ----- ----- ----- ----- ----- -----
            %   * = HW region
            %   Remarks:    1 < ymin < ymax < ny (GV.children.cartDims(2))
            %               1 < zmin < zmax < nz (GV.layers.num)
            %   regionIndices: [ymin, ymax, zmin, zmax]
            %
            hw.GVOI           = G;
            hw.well           = well;
            hw.regionIndices  = regionIndices;
            hw.checkRegionIndices();
        end
        
        function checkRegionIndices(hw)
            % Check whether the region indices exceed the Cartesian reigon
            % dimension of the VOI grid
            [nx, ny, nz] = hw.assignCartDimsOfVOIGrid();
            assert(hw.well.segmentNum == nx);
            [ymin, ymax, zmin, zmax] = hw.assignRegionIndices();
            assert(1 < ymin & ymin < ymax & ymax < ny, ...
                '1 < Indices(1) < Indices(2) < %2.0f  is not satisfied', ny)
            assert(1 < zmin & zmin < zmax & zmax < nz, ...
                '1 < Indices(3) < Indices(4) < %2.0f  is not satisfied', nz)
        end
        
        function packed = allInfoOfRegion(hw, varargin)
            % Get all information of the region, consisting of cells,
            % layer-faces, nodes of layer-faces, boundary nodes, indices of
            % vertices.
            packed.cells = hw.cellsOfRegion();
            packed.nodes = hw.nodesOfRegion();
            packed.faces = hw.facesOfRegion(packed.cells, packed.nodes);
            packed.bdyNodes  = hw.bdyNodesOfRegion();
            packed.vertexID  = hw.IDOfFourVertices();
        end
        
        function c = cellsOfRegion(hw, varargin)
            % Get cells of HW region in VOI grid
            G  = hw.GVOI;
            [nx, ny, nz] = hw.assignCartDimsOfVOIGrid();
            [ymin, ymax, zmin, zmax] = hw.assignRegionIndices();
            c_y = (ymin : ymax)';
            c_yz = repmat(c_y, 1, nz);
            c_yz = bsxfun(@plus, c_yz, (0:nz-1) * G.cells.num / nz);
            c_layers = (zmin : zmax)';
            c_yz = c_yz(:, c_layers);
            c_yz = c_yz(:);
            c = arrayfunUniOut(@(x) c_yz + x * ny, (0:nx-1)');
        end
        
        function n = nodesOfRegion(hw, varargin)
            % Get nodes of layer-faces of HW region in VOI grid
            [nx, ny] = hw.assignCartDimsOfVOIGrid();
            n_yz = hw.getNodesSingleSurface();
            n = arrayfunUniOut(@(x)n_yz(:) + x * (ny+1), (0:nx)');
        end
        
        function f = facesOfRegion(hw, c, n, varargin)
            % Get layer-faces of HW region in VOI grid
            G = hw.GVOI;
            faceFun = @(c)gridCellFaces(G, c);
            nodeFun = @(f)gridFaceNodes(G, f);
            f = cell(length(n), 1);
            c = [c; c(end)];
            for k = 1 : numel(f)
                ck = c{k};
                fk = arrayfunUniOut(faceFun, ck);
                f{k} = zeros(numel(ck),1);
                for j = 1 : length(ck)
                    fkj = fk{j};
                    nds = arrayfunUniOut(nodeFun, fkj);
                    idx = cellfun(@(x)all(ismember(x, n{k})), nds);
                    f{k}(j) = fkj(idx);
                end
            end
        end
        
        function bn = bdyNodesOfRegion(hw, varargin)
            % Get boundary nodes of HW region in VOI grid
            [nx, ny] = hw.assignCartDimsOfVOIGrid();
            n_yz = hw.getNodesSingleSurface(hw);
            bn_yz = [n_yz(1, :)'; n_yz(2:end-1, end); ...
                n_yz(end, end:-1:1)'; n_yz(end-1:-1:2, 1)];
            bn = arrayfunUniOut(@(x) bn_yz + x * (ny+1), (0:nx)');
        end
        
        function vxID = IDOfFourVertices(hw)
            % Get the indices of four vertices in boundary nodes
            n_yz = hw.getNodesSingleSurface();
            vx    = [n_yz(1,1); n_yz(1,end); n_yz(end,end); n_yz(end,1)];
            bn_yz = [n_yz(1, :)'; n_yz(2:end-1, end); ...
                n_yz(end, end:-1:1)'; n_yz(end-1:-1:2, 1)];
            vxID = find(ismember(bn_yz, vx));
        end
        
        function c = cartCellsOfVOIGrid(hw, varargin)
            % Get cells of VOI grid in Cartesian region
            G  = hw.GVOI;
            [nx, ny, nz] = hw.assignCartDimsOfVOIGrid();
            c = bsxfun(@plus, (1:ny), G.cells.num/nz * (0:nz-1)');
            c = arrayfunUniOut(@(i)c(:) + (i-1)*ny, (1:nx)');
        end
        
        function n = getNodesSingleSurface(hw, varargin)
            % Get nodes of HW region in VOI grid on single surface
            G  = hw.GVOI;
            [~, ~, nz] = hw.assignCartDimsOfVOIGrid();
            [ymin, ymax, zmin, zmax] = hw.assignRegionIndices();
            n_y  = (ymin : ymax+1)';
            n_yz = repmat(n_y, 1, nz+1);
            n_yz = bsxfun(@plus, n_yz, (0:nz) * G.nodes.num / (nz + 1));
            n_layers = (zmin : zmax+1)';
            n_yz = n_yz(:, n_layers);
            n  = n_yz';
        end
        
        function [ymin, ymax, zmin, zmax] = assignRegionIndices(hw)
            % Assign the region indices
            Ind = hw.regionIndices;
            [ymin, ymax, zmin, zmax] = deal(Ind(1), Ind(2), Ind(3), Ind(4));
        end
        
        function [nx, ny, nz] = assignCartDimsOfVOIGrid(hw)
            % Assign Cartesian dimensions of VOI grid
            G  = hw.GVOI;
            G2 = G.surfGrid;
            [nx, ny, nz] = deal(G2.cartDims(1),G2.cartDims(2),G.layers.num);
        end
        
        function GW = ReConstructToRadialGrid(hw, radPara)
            % Reconstruct VOI grid in HW region to layered radial grid.
            % Two types of grid lines are provided:
            % 'pureCircular' : The radial grid lines are pure circular
            % 'gradual'      : The radial grid lines vary from the circular
            %                  line to the rectangular line of a specified
            %                  box gradually
            % PARAMETERS:
            %  GV      - The layered VOI grid, built by 
            %            'VolumeOfInterest.ReConstructGrid'
            %  radPara - Parameters for generating the radial grid
            %            The type 'pureCircular' requires following fields:
            %               'maxRadius': Max radius of the radial grid
            %               'nRadCells': Number of radial cells
            %            The type 'gradual' requires following fields:
            %               'boxRatio' : Size ratio of the rectangular box 
            %                            to the outer boundary, 2x1 double, 
            %                            [yRatio, zRatio]
            %               'nRadCells': Number of radial cells, 2x1 double,    
            %                            [inbox, outbox]
            %               'pDMult'   : Multiplier of pD of the outer-most 
            %                            angular line to PD of wellbore
            %                            line, the larger the pDMult is,
            %                            the outer-most line is closer to
            %                            the box
            %              'offCenter':  Whether the well is off-center in 
            %                            the radial grid
            %                            
            % RETURNS:
            %  GW - Layered radial HW grid
            
            fprintf(' -- Reconstructing the VOI grid to radial HW grid\n')
            
            % Geometrical info of HW reigon
            packed = hw.allInfoOfRegion();
            [pSurfs, pSurfXY, wellbores] = ...
                generateHWGridNodes(hw.GVOI, packed, hw.well, radPara);
            
            % Construct the 2D grid
            p   = pSurfXY{1};
            nA  = unique(cellfun(@numel, packed.bdyNodes));
            assert(numel(nA)==1);
            [gW, t] = buildRadialGrid(p, nA, sum(radPara.nRadCells));
            gW.nodes.boundary = (gW.nodes.num-nA+1 : gW.nodes.num)';
            
            % Rewrite the radial dimensions, in order to be compatible with
            % 'computeRadTransFactor'
            switch radPara.gridType
                case 'pureCircular'
                    % The outer-most radial cells are not 'real radial cells'
                    gW.radDims = [nA, radPara.nRadCells-1, 1];
                case 'gradual'
                    % The radial cells outside the box are not 'real radial
                    % cells'
                    gW.radDims = [nA, radPara.nRadCells];
                otherwise
                    error([mfilename, ': Unknown radial grid type'])
            end
            
            % Extrude the 2D grid to 3D grid
            GW = makeLayeredGridNWM(gW, pSurfs, 'connectivity', t);
            GW.radDims = [gW.radDims, GW.layers.num];
            GW.layers.refinement = ones(GW.layers.num, 1);
            GW.wellbores = wellbores;
            GW.layers.coordsXY = pSurfXY;
            GW.parentInfo = packed;
        end
        
        function showWellRegionInVOIGrid(hw, varargin)
            % Show the well region in the VOI grid
            G   = hw.GVOI;
            pW  = hw.well.trajectory;
            opt = struct('showWellRgionCells', true);
            opt = merge_options(opt, varargin{:});
            cCart = hw.cartCellsOfVOIGrid();
            cReg  = hw.cellsOfRegion();
            cRes = cCart{1}(~ismember(cCart{1}, cReg{1}));
            figure, hold on
            plotGrid(G, cRes, 'facecol', rand(3,1))
            if opt.showWellRgionCells
                plotGrid(G, cReg{1}, 'facecol', rand(3,1))
                plot3(pW(1,1), pW(1,2), pW(1,3), ...
                    'rs', 'markerfacecol', 'r', 'markersize', 10)
            else
            end
            axis equal tight
        end

        function plotRegionCells(hw, packed)
            % Plot cells inside the HW region
            figure, hold on, axis off
            arrayfun(@(x)plotGrid(hw.GVOI, packed.cells{x}, 'facecolor',...
                rand(3,1)), (1:numel(packed.cells))')
            view(3)
            title('Cells inside the region')
        end
        
        function plotRegionLayerFaces(hw, packed)
            % Plot layer-faces of the HW region
            figure, hold on, axis off
            arrayfun(@(x)plotFaces(hw.GVOI, packed.faces{x}, 'facecolor',...
                rand(3,1)), (1:numel(packed.faces))')
            view(3)
            title('Layer-faces of the region')
        end
    end
end