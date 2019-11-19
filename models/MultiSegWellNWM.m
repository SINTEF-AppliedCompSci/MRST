classdef MultiSegWellNWM < NearWellboreModel
% Derived class for collecting simulation data of the near-wellbore model
% coupling with multi-segment well model

    properties
       wellboreGrid 
    end
    
    methods
        function msw = MultiSegWellNWM(subGrids, deck, well, varargin)
            msw = msw@ NearWellboreModel(subGrids, deck, well, varargin);
            assert( isfield(well, 'isMS') && well.isMS, ...
                'The input well is not a multi-segment well' )
            msw.wellboreGrid = msw.buildWellboreGrid();
        end
        
        function model = setupSimModel(nwm, rock, T_all, N_all)
            % Setup simulation model passed to ad-blackoil solver for the
            % global grid (the multi-segment well model now only support
            % 'ThreePhaseBlackOilModel')
            % rock:  Rock of global grid
            % T_all: Full transmissibility
            % N_all: Neighbor ship
            gravity reset on
            G    = nwm.gloGrid;
            f    = nwm.fluid;
            % Internal connections
            intCon = all(N_all, 2);
            N = N_all(intCon, :);
            T = T_all(intCon);
            % Phase component
            ph = nwm.getPhaseFromDeck();
            assert(all(ph.wat & ph.oil & ph.gas), ['The multi-segment '...
                'well model now only support ''ThreePhaseBlackOilModel'''])
            model = ThreePhaseBlackOilModel(G, rock, f, 'water', ph.wat, ...
                'oil', ph.oil, 'gas', ph.gas, 'vapoil', ph.vapo, 'disgas', ph.disg);
            % Reset the operators
            model.operators = setupOperatorsTPFA(G, rock, ...
                'neighbors', N, 'trans', T);
            model.operators.T_all = T_all;
        end

        function schedule = getSimSchedule(msw, model, varargin)
            % Get the multi-segment well simulation schedule from deck and
            % node/segment definition
            schedule = getSimSchedule@NearWellboreModel(msw, model);
            nodes = msw.generateNodes();
            segs = msw.generateSegments();
            for i = 1 : numel(schedule.control)
                W0 = schedule.control(i).W;
                name = vertcat(W0.name);
                ii = all(name == msw.well.name, 2);
                WOthers = W0(~ii);
                W = W0(ii);
                W = convert2MSWell(W, ...
                    'cell2node'   , nodes.cell2node, ...
                    'connDZ'      , [], ...
                    'nodeDepth'   , nodes.depth, ...
                    'topo'        , segs.topo, ...
                    'segLength'   , segs.length, ...
                    'segRoughness', [], ...
                    'segFlowModel', [], ...
                    'segType'     , [], ...
                    'segDiam'     , segs.diam, ...
                    'G'           , [], ...
                    'vol'         , nodes.vol);
                W.segments.roughness = segs.roughness;
                W.segments.flowModel = segs.flowModel;
                % Redefine the reference depth
                dz = W.refDepth - nodes.depth(1);
                W.refDepth = nodes.depth(1);
                W.dZ = W.dZ + dz;
                WNew = combineMSwithRegularWells(WOthers, W);
                schedule.control(i).W = WNew;
            end
        end
        
        function gW = buildWellboreGrid(msw)
            % Build grid for the space inside wellbore
            [~, ~, GW] = msw.assignSubGrds();
            % Get the wall and screen nodes
            wellbores = GW.wellbores;
            pW = cellfunUniOut(@(x)x.wall.coords, wellbores);
            nA = GW.radDims(1);
            assert(all( cellfun(@(x)size(x,1)==nA, pW) ));
            % Build the radial grid
            p = pW{1}(:, [1,2]);
            t = {(1:nA)'};
            % Build the wellbore grid
            g = tessellationGrid(p, t);
            g.nodes.coords = pW{1};
            gW = makeLayeredGridNWM(g, pW, 'connectivity', t);
            gW.radDims = [nA, 1, gW.layers.num];
        end

        function nodes = generateNodes(msw)
            % Generate nodes information from the wellbore grid for 
            % multi-segment well
            gW = msw.wellboreGrid;
            nA = gW.radDims(1);
            n = msw.well.openedSegs;
            n = repmat(n, nA, 1);
            n = n(:);
            wc = msw.getWellCells();
            cell2node = sparse(n, (1:numel(wc))', 1, gW.cells.num, numel(wc));
            nodes = struct(...
                'coords'    , gW.cells.centroids, ...
                'depth'     , gW.cells.centroids(:,3), ...
                'vol'       , gW.cells.volumes, ...
                'dist'      , nan(gW.cells.num,1), ...
                'cell2node' , cell2node, ...
                'resCells'  , wc);
        end
        
        function segs = generateSegments(msw)
            % Generate segment information from the wellbore grid for
            % multi-segment well
            gW = msw.wellboreGrid;
            % Topology
            f   = all(gW.faces.neighbors,2);
            t   = gW.faces.neighbors(f,:);
            t   = sort(t, 2);
            % Length
            dxyz = gW.cells.centroids(t(:,1), :) - gW.cells.centroids(t(:,2), :);
            L  = sqrt( sum(dxyz.^2, 2) );
            % Roughness
            roughness = msw.well.roughness;
            roughness = (roughness(1:end-1)+roughness(2:end))/2;
            % Diameter
            D  = 2*msw.well.radius(2:end-1);
            % Area
            area = gW.faces.areas(f);
            % Segment number
            ns = numel(L);
            fm = @(v, rho, mu)wellBoreFriction(v, rho, mu, D, L, roughness, 'massRate');
            segs = struct(...
                'length'   , L, ...
                'roughness', roughness, ...
                'diam'     , D, ...
                'topo'     , t, ...
                'area'     , area, ...
                'num'      , ns, ...
                'flowModel', fm);
        end
        
        function plotCell2Node(msw, nodes, N)
            % Plot reservoir cells connected to node N
            G = msw.gloGrid;
            coords = nodes.coords;
            figure, hold on
            c2n = nodes.cell2node;
            wc = nodes.resCells;
            if size(N, 1) > size(N,2); N = N'; end;
            for i = N
                idx = find(c2n(i, :));
                for j = 1 : numel(idx)
                    p = [coords(i,:); G.cells.centroids(wc(idx(j)), :)];
                    plot3(p(:,1), p(:,2), p(:,3), '.-', 'color', [1, 0, 0])
                end
                plotGrid(G, wc(idx), 'facecolor', 'none')
            end
        end
        
    end
end