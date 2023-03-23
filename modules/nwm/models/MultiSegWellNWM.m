classdef MultiSegWellNWM < NearWellboreModel
% Derived class for generating necessary data structures passed to the mrst
% AD simulators for the hybrid grid of near-wellbore model coupling with 
% the multi-segment well model

    properties
       wellboreGrid  % 1D 'wellbore grid' in the void wellbore space which
                     % conforms with the reservoir grid
    end
    
    methods
        function msw = MultiSegWellNWM(subGrids, deck, well, varargin)
            msw = msw@NearWellboreModel(subGrids, deck, well, varargin);
            assert( isfield(well, 'isMS') && well.isMS, ...
                'The input well is not a multi-segment well' )
            msw.wellboreGrid = msw.buildWellboreGrid();
        end
        
        function model = setupSimModel(nwm, rock, T_all, N_all)
            % Setup simulation model passed to ad-blackoil simulator for 
            % the global grid (the multi-segment well model now only
            % supports the 'ThreePhaseBlackOilModel')
            % rock:  Rock of the global grid
            % T_all: Full transmissibility
            % N_all: Neighborship of all connections
            gravity reset on
            G    = nwm.gloGrid;
            f    = nwm.fluid;
            % Internal connections
            intCon = all(N_all, 2);
            N = N_all(intCon, :);
            T = T_all(intCon);
            % Phase components
            ph = nwm.getPhaseFromDeck();
            assert(all(ph.wat & ph.oil & ph.gas), ['The multi-segment '...
                'well model now only support ''ThreePhaseBlackOilModel'''])
            model = ThreePhaseBlackOilModel(G, rock, f, 'water', ph.wat, ...
                'oil', ph.oil, 'gas', ph.gas, 'vapoil', ph.vapo, 'disgas', ph.disg);
            % Reset the operators
            model.operators = setupOperatorsTPFA(G, rock, 'neighbors', N, 'trans', T);
            model.operators.N_all = N_all;
            model.operators.T_all = T_all;
            
            % Aquifer model
            [hasAquifer, output] = nwm.handleAquifers();
            if hasAquifer
                model.AquiferModel = AquiferModel(model, ...
                    output.aquifers, output.aquind, output.aquiferprops, output.initval);
            end
        end

        function schedule = getSimSchedule(msw, model, varargin)
            % Get the multi-segment well simulation schedule from deck and
            % node/segment definitions
            opt = struct('returnMS', true, 'refDepthFrom', 'topNode');
            opt = merge_options(opt, varargin{:});
            schedule = getSimSchedule@NearWellboreModel(msw, model, ...
                'refDepthFrom', 'topNode');
            if ~opt.returnMS
                return
            end
            nodes = msw.generateNodes();
            segs = msw.generateSegments();
            for i = 1 : numel(schedule.control)
                W0 = schedule.control(i).W;
                ii = arrayfun(@(w)strcmp(msw.well.name, w.name), W0);
                WRegular = W0(~ii);
                W = W0(ii);
                W = convert2MSWell(W, ...
                    'cell2node'   , nodes.cell2node, ...
                    'connDZ'      , W.dZ, ...
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
                if numel(WRegular) > 0
                    WNew = combineMSwithRegularWells(WRegular, W);
                else
                    WNew = W;
                end
                schedule.control(i).W = WNew;
            end
        end
        
        function gW = buildWellboreGrid(msw)
            % Build grid for the void space inside wellbore
            [~, ~, GW] = msw.assignSubGrds();
            % Get the borewall (casing) nodes and connectivity list
            wellbores = GW.wellbores;
            pW = cellfunUniOut(@(x)x.wall.coords, wellbores);
            nA = GW.radDims(1);
            assert(all( cellfun(@(x)size(x,1)==nA, pW) ));
            p = pW{1}(:, [1,2]);
            t = {(1:nA)'};
            % Build the wellbore grid
            g = tessellationGrid(p, t);
            g.nodes.coords = pW{1};
            gW = makeLayeredGridNWM(g, pW, 'connectivity', t);
            gW.radDims = [nA, 1, gW.layers.num];
        end

        function nodes = generateNodes(msw)
            % Generate node definitions from the wellbore grid for 
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
            % Generate segment definitions from the wellbore grid for
            % multi-segment well
            gW = msw.wellboreGrid;
            % Topology
            f   = find( all(gW.faces.neighbors,2) );
            t   = gW.faces.neighbors(f,:);
            t   = sort(t, 2);
            % Length
            dxyz1 = gW.cells.centroids(t(:,1), :) - gW.faces.centroids(f, :);
            L1    = sqrt( sum(dxyz1.^2, 2) );
            dxyz2 = gW.cells.centroids(t(:,2), :) - gW.faces.centroids(f, :);
            L2    = sqrt( sum(dxyz2.^2, 2) );
            L     = L1 + L2;
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
        
        function plotSegments(msw, nodes, segs, S)
            % Plot the nodes and reservoir cells assoicated with segment S
            G = msw.gloGrid;
            coords = nodes.coords;
            topo = segs.topo(S, :);
            figure, hold on
            c2n = nodes.cell2node;
            wc = nodes.resCells;
            N = unique(topo(:));
            for i = N'
                idx = find(c2n(i, :));
                for j = 1 : numel(idx)
                    p = [coords(i,:); G.cells.centroids(wc(idx(j)), :)];
                    plot3(p(:,1), p(:,2), p(:,3), 's-', 'color', ...
                        [1, 0, 0], 'linewidth', 1)
                end
                plotGrid(G, wc(idx), 'facecolor', 'none')
            end
            for k = 1 : size(topo,1)
                p = coords(topo(k, :), :);
                plot3(p(:,1), p(:,2), p(:,3), 's-', 'color', ...
                    [0, 0, 1], 'linewidth', 1)
            end
        end
        
    end
end