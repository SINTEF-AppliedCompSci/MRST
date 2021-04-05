classdef NearWellboreModel
% Class for generating necessary data structures passed to the mrst AD 
% simulators for the hybrid grid of near-wellbore model

    properties
        subGrids    % Subgrids {Corner-point grid (CPG), VOI grid, HW grid}
        inputDeck   % ECLIPSE-style input deck of CPG
        well        % Structure of well information
        gloGrid     % The global hybrid grid
        fluid       % AD-solver fluid from ECLIPSE-style input deck
    end

    methods
        function nwm = NearWellboreModel(subGrids, deck, well, varargin)
            checkDeck(deck);
            nwm.subGrids  = subGrids;
            nwm.inputDeck = deck;
            nwm.well      = well;
            nwm.gloGrid   = nwm.validateGlobalGrid();
            nwm.fluid     = nwm.setupFluid();
        end

        function G = validateGlobalGrid(nwm, varargin)
            % Validate the global hybrid grid by subgrids
            % Get updated subgrids (cells in VOI and HW region are removed)
            GCu = nwm.updateCPG();
            GVu = nwm.updateVOIGrid();
            GWu = nwm.updateHWGrid();
            % Combine the updated subgrids
            G = assembleGrids({GCu, GVu, GWu});
        end

        function f = setupFluid(nwm)
            % Initialize AD fluid from input deck
            deck  = nwm.inputDeck;
            f = initDeckADIFluid(deck);
        end
        
        function varargout = assignInputSubGrds(nwm, varargin)
            % Assign subgrids from input
            assert(nargout <= numel(nwm.subGrids))
            if nargout == 1
                varargout = {nwm.subGrids};
            else
                varargout = nwm.subGrids(1:nargout);
            end
        end

        function varargout = assignSubGrds(nwm, varargin)
            % Assign updated subgrids from the global grid
            G = nwm.gloGrid;
            assert(nargout <= numel(G.subGrids))
            if nargout == 1
                varargout = {G.subGrids};
            else
                varargout = G.subGrids(1:nargout);
            end
        end

        function varargout = assignInputSubGrdTypes(nwm, varargin)
            % Assign types of input subgrids
            assert(nargout <= numel(nwm.subGrids))
            types = {'CPG', 'VOI Grid', 'HW Grid'};
            if nargout == 1
                varargout = {types};
            else
                varargout = types(1:nargout);
            end
        end

        function varargout = assignSubGrdTypes(nwm, varargin)
            % Assign types of updated subgrids
            G = nwm.gloGrid;
            assert(nargout <= numel(G.subGrids))
            types = {'Updated CPG', 'Updated VOI Grid', 'Updated HW Grid'};
            if nargout == 1
                varargout = {types};
            else
                varargout = types(1:nargout);
            end
        end

        function [G, rock, f, model, schedule, initState] = packedSimData(nwm, rockW, varargin)
            % Obtain all necessary simulation data structures of the hybrid
            % grid in near-wellbore model
            % Global grid
            G = nwm.gloGrid;
            % AD fluid
            f = nwm.fluid;
            % Rocks
            rockC = nwm.getCPGRockFromDeck();
            rockV = nwm.getVOIRocksByInterp();
            rock  = nwm.getGlobalRock({rockC, rockV, rockW});
            % Simulation model
            T = nwm.getTransGloGrid(rock);
            intXn = nwm.computeIntxnRelation();
            nnc = nwm.generateNonNeighborConn(intXn, rock, T);
            G.nnc = nnc;
            [T_all, N_all] = nwm.assembleTransNeighbors(T, nnc);
            model = nwm.setupSimModel(rock, T_all, N_all);
            % Schedule
            schedule = nwm.getSimSchedule(model, varargin{:});
            % Initial state
            initState = nwm.getInitState(model);
        end
        
        function [GC, rockC, f, modelC, scheduleC, initStateC] = packedCPGSimData(nwm)
            % Obtain all necessary simulation data structures of the CPG
            % CPG
            [GC, ~] = nwm.assignInputSubGrds();
            % AD fluid
            f = nwm.fluid;
            % Rocks
            rockC = nwm.getCPGRockFromDeck();
            % Simulation model
            modelC = nwm.setupCPGSimModel();
            % Schedule
            scheduleC = nwm.getCPGSimSchedule(modelC);
            % Init state
            initStateC = nwm.getInitState(modelC);
        end
        
        function rock = getGlobalRock(nwm, rocks)
            % Get the rock for the global grid
            % rocks = {rockC, rockV, rockW}
            % -------------------------------------------------------------
            % | Rock  | Grid  | Source        | Permeability | Anisotropy |
            % |       |       |               | coordinate   |            |
            % |-----------------------------------------------------------|
            % | rockC | Input | Input deck    | Local        | Yes        |
            % |       | GC    |               |              |            |
            % |-----------------------------------------------------------|
            % | rockV | Input | Interpolation | Global       | Yes        |
            % |       | GV    | of rockC      |              |            |
            % |-----------------------------------------------------------|
            % | rockW | Input | User-defined  | Global       | No         |
            % |       | GW    |               |              |            |
            % -------------------------------------------------------------
            for j = 1 : numel(rocks)
                if ~isfield(rocks{j}, 'ntg')
                    rocks{j}.ntg = ones(size(rocks{j}.poro));
                end
            end
            % Map the rocks from input subgrids to global grid
            G = nwm.gloGrid;
            mapc = nwm.cellMapFromInputSubGrdsToGloGrd();
            fn = {'perm', 'poro', 'ntg'};
            for i = 1 : numel(fn)
                for j = 1 : numel(rocks)
                    rock.(fn{i})(mapc{j}(:,2), :) = ...
                        rocks{j}.(fn{i})(mapc{j}(:,1), :);
                end
                assert(size(rock.(fn{i}), 1) == G.cells.num)
            end
        end
        
        function rockC = getCPGRockFromDeck(nwm)
            % Get the rock of input CPG from input deck
            deck = nwm.inputDeck;
            [GC, ~] = nwm.assignInputSubGrds();
            rockC  = initEclipseRock(deck);
            rockC  = compressRock(rockC, GC.cells.indexMap);
        end

        function rockV = getVOIRocksByInterp(nwm, varargin)
            % Get the rock of input VOI grid by interpolation of CPG rock
            % Optional:
            %  'InterpMethod', same with opitions in 'griddata':
            %  'linear' (default) | 'nearest' | 'natural' | 'cubic' | 'v4'
            opt = struct('InterpMethod', 'linear');
            opt = merge_options(opt, varargin{:});

            [GC, GV] = nwm.assignInputSubGrds();
            rockC    = nwm.getCPGRockFromDeck();
            grdecl   = nwm.getGrdEclFromDeck();
            
            % Layer indices to determine the corresponding layers in packed 
            % data of the VOI grid layers
            layerID = rldecode((1:length(GV.layers.refinement))', ...
                GV.layers.refinement);
            assert(length(layerID) == GV.layers.num)

            % Get CPG cell centers and cellFace centers
            if isfield(grdecl, 'COORD')
                [cCenters, cfCenters] = computeCpGeometry(GC, grdecl);
            else
                cCenters  = GC.cells.centroids;
                cfCenters = GC.faces.centroids(GC.cells.faces(:,1), :);
            end
            
            % Get the properties by interpolation
            rockV.perm = zeros(GV.cells.num,3);
            rockV.poro = zeros(GV.cells.num,1);
            rockV.ntg  = zeros(GV.cells.num,1);
            packed = GV.parentInfo;
            for layer =  1 : GV.layers.num
                cellsC = packed.boxCells{layerID(layer)};
                % Perm, loc --> glo -----------------------------------
                [ux, uy, uz] = getUnitDisVectors(GC, cfCenters, cellsC);
                perm_loc = rockC.perm(cellsC,:);
                perm_glo = bsxfun(@times, perm_loc(:,1), ux) + ...
                    bsxfun(@times, perm_loc(:,2), uy) + ...
                    bsxfun(@times, perm_loc(:,3), uz);
                xx = cCenters(cellsC,1);
                yy = cCenters(cellsC,2);
                cellsV = find( GV.cells.layers == layer );
                xq = GV.cells.centroids(cellsV,1);
                yq = GV.cells.centroids(cellsV,2);
                permV = zeros(length(cellsV), 3);
                for dim = 1 : 3
                    permV(:,dim) = griddata(xx, yy, perm_glo(:,dim), ...
                        xq, yq, opt.InterpMethod);
                end
                rockV.perm(cellsV, :) = permV;

                % poro ----------------------------------
                poroC = rockC.poro(cellsC);
                poroV = griddata(xx, yy, poroC, xq, yq, opt.InterpMethod);
                rockV.poro(cellsV) = poroV;

                % ntg -----------------------------------
                ntgC = rockC.ntg(cellsC);
                ntgV = griddata(xx, yy, ntgC, xq, yq, opt.InterpMethod);
                rockV.ntg(cellsV) = ntgV;
            end
        end

        function varargout = assignSubRocks(nwm, rock)
            % Assign rocks from the global rock for updated subgrids 
            mapc = nwm.cellMapFromSubGrdsToGloGrd();
            fn = fieldnames(rock);
            subRocks = cell(numel(mapc),1);
            for i = 1 : numel(subRocks)
                for j = 1 : numel(fn)
                    subRocks{i}.(fn{j})(mapc{i}(:, 1), :) = ...
                        rock.(fn{j})(mapc{i}(:, 2), :);
                end
            end
            assert(nargout <= numel(nwm.subGrids))
            if nargout == 1
                varargout = {subRocks};
            else
                varargout = subRocks(1:nargout);
            end
        end

        function [T_all, N_all] = assembleTransNeighbors(nwm, T, nnc)
            % Assemble transmissibility and neighborship for the simulation
            % model
            G = nwm.gloGrid;
            T_all = [T;                 nnc.T];
            N_all = [G.faces.neighbors; nnc.cells];
            assert(numel(T_all) == size(N_all,1))
        end
        
        function T = getTransGloGrid(nwm, rock)
            % Compute the transmissibility for the global grid
            % rock: rock of the global grid
            % The transmissibility consists of: 
            % Transmissibility of updated [CPG, VOI grid, and HW grid]
            
            % Half transmissibility of updated CPG
            hTC = nwm.computeCPGHalfTrans(rock);
            % Half transmissibility of updated VOI grid
            hTV = nwm.computeVOIGrdHalfTrans(rock);
            % Half transmissibility of updated HW grid
            hTW = nwm.computeHWGrdHalfTrans(rock);
            % Combine the half transmissibilities
            hT = [hTC; hTV; hTW];
            % Get full transmissibility, corresponding to G.faces.neighbors
            G = nwm.gloGrid;
            cf = G.cells.faces(:,1);
            assert( numel(hT) == numel(cf) )
            nf = G.faces.num;
            T  = 1 ./ accumarray(cf, 1./hT, [nf, 1]);
        end
        
        function hT = computeCPGHalfTrans(nwm, rock)
            % Compute half transmissibility of the updated CPG
            % -----------------------------------------------
            % | Flow          | Permeability   | Anisotropy |
            % | approximation | coordinate     |            |
            % |---------------------------------------------|
            % | Linear        | Local          | Yes        |
            % -----------------------------------------------
            
            [GCu, ~] = nwm.assignSubGrds();
            [rockCu, ~] = nwm.assignSubRocks(rock);
            grdecl = nwm.getGrdEclFromDeck();
            if isfield(grdecl, 'COORD')
                [cCenters, cFCenters] = computeCpGeometry(GCu, grdecl);
                hT = computeTrans(GCu, rockCu, 'K_system', 'loc_xyz', ...
                    'cellCenters', cCenters, ...
                    'cellFaceCenters', cFCenters);
            elseif isfield(grdecl, 'DX')
                hT = computeTrans(GCu, rockCu);
            else
                error([mfilename, ': Unknown deck grid input'])
            end
        end
        
        function hT = computeVOIGrdHalfTrans(nwm, rock)
            % Compute half transmissibility of the updated VOI Grid
            % -----------------------------------------------
            % | Flow          | Permeability   | Anisotropy |
            % | approximation | coordinate     |            |
            % |---------------------------------------------|
            % | Linear        | Global         | No         |
            % -----------------------------------------------
            [~, GVu] = nwm.assignSubGrds();
            [~, rockVu] = nwm.assignSubRocks(rock);
            hT = computeTrans(GVu, rockVu);
        end

        function hT = computeHWGrdHalfTrans(nwm, rock)
            % Compute half transmissibility of the updated HW Grid
            % -----------------------------------------------
            % | Flow          | Permeability   | Anisotropy |
            % | approximation | coordinate     |            |
            % |---------------------------------------------|
            % | Radial        | Global         | Yes        |
            % -----------------------------------------------
            [~, ~, GWu] = nwm.assignSubGrds();
            [~, ~, rockWu] = nwm.assignSubRocks(rock);
            % Compute the linear transmissibility first
            hT = computeTrans(GWu, rockWu);
            % Get the radial transmissibility factor
            ft = nwm.getRadTransFactors();
            assert(numel(ft) == numel(GWu.cells.faces(:,1)))
            % Compute the radial transmissibility
            DZ = arrayfun(@(c)getDZ(GWu, c), (1:GWu.cells.num)');
            DZ = rldecode(DZ, diff(GWu.cells.facePos));
            perm = rockWu.perm(:,1);
            perm = rldecode(perm, diff(GWu.cells.facePos));
            hT_rad = perm .* DZ .* ft;
            % Assgin the radial transmissibility
            isRad = ~isnan(hT_rad);
            hT(isRad) = hT_rad(isRad);
        end
        
        function ft = getRadTransFactors(nwm)
            % Get the radial half transmissibility factors for the HW grid
            % Get the factors of the 2D surface grids first
            fprintf(' -- Computing the radial transmissibility factors\n')
            [~, ~, GWu] = nwm.assignSubGrds();
            gW    = GWu.surfGrid;
            pXYs  = GWu.layers.coordsXY;
            % Skin factors of segments
            s_seg = nwm.well.skinFactor;
            if size(s_seg,1) > size(s_seg,2); s_seg = s_seg'; end
            % Assign to all surface grids
            skins = [s_seg(1), (s_seg(1:end-1)+ s_seg(2:end))/2, s_seg(end)];
            assert(numel(s_seg) == nwm.well.segmentNum, ...
                'The skin factors should be given for all segments')
            pW = [0, 0];
            ft = arrayfunUniOut(@(i)computeRadTransFactor(gW, pW, skins(i), ...
                'nodeCoords', pXYs{i}), (1:numel(pXYs))');
            % Extend the factor to layered HW grid
            ft = cellfun(@(x,y)(x+y)/2, ft(1:end-1), ft(2:end), ...
                'UniformOutput', false);
            for k = 1 : length(ft)
                tmp = reshape(ft{k}, 4, []);
                tmp = [tmp; nan(2, size(tmp,2))];
                tmp = tmp(:);
                ft{k} = tmp;
            end
            ft = cell2mat(ft);
        end
        
        function nnc = generateNonNeighborConn(nwm, intXn, rock, T)
            % Generate the non-neighbor connections (NNCs)
            % intXn: Boundary intersection relations
            % rock:  Rock of global grid
            % T:     Fully transmissibility of global grid
            nnc_nmf = nwm.nncOfNonMatchingBoundaries(intXn, rock, T);
            nnc_mf  = nwm.nncOfMatchingBoundaries(intXn, rock, T);
            nnc = struct('cells', [], 'T', [], 'hT', []);
            fn = fieldnames(nnc);
            for i = 1 : numel(fn)
                nnc.(fn{i}) = vertcat(nnc_nmf.(fn{i}), nnc_mf.(fn{i}));
            end
        end
        
        function nnc = nncOfNonMatchingBoundaries(nwm, intXn, rock, T)
            % Generate non-neighbor connections (NNCs) arised from the 
            % non-matching boundaries
            % intXn: Boundary intersection relations
            % rock:  Rock of global grid
            % T:     Fully transmissibility of global grid
            G  = nwm.gloGrid;
            nmf = intXn.nonMatchingFaces;
            % Centers of subfaces
            fc = nmf(:, 4:6);
            % Area normals of subfaces
            N  = nmf(:, 7:9);
            num = size(nmf,1);
            nnc = struct('cells', zeros(num, 2), ...
                'T', zeros(num, 1), 'hT', zeros(num, 2));
            for j = 1 : 2
                faces = nmf(:,j);
                assert(all( ~all(G.faces.neighbors(faces, :), 2) ))
                cells = sum(G.faces.neighbors(faces, :), 2);
                K = rock.perm(cells, :);
                D = fc - G.cells.centroids(cells,:);
                hT = K(:,1).*D(:,1).*N(:,1) + K(:,2).*D(:,2).*N(:,2) + ...
                    K(:,3).*D(:,3).*N(:,3);
                hT = abs(hT) ./ sum(D .* D, 2);
                nnc.cells(:, j) = cells;
                nnc.hT(:, j)    = hT;
            end
            nnc.T  = 1./(1./nnc.hT(:,1) + 1./nnc.hT(:,2));
        end
        
        function nnc = nncOfMatchingBoundaries(nwm, intXn, rock, T)
            % Generate non-neighbor connections (NNCs) arised from the 
            % matching boundaries
            % intXn: Boundary intersection relations
            % rock:  Rock of global grid
            % T:     Fully transmissibility of global grid
            G  = nwm.gloGrid;
            mf = intXn.matchingFaces;
            num = size(mf,1);
            nnc = struct('cells', zeros(num, 2), ...
                'T', zeros(num, 1), 'hT', zeros(num, 2));
            for j = 1 : 2
                faces = mf(:,j);
                assert(all( ~all(G.faces.neighbors(faces, :), 2) ))
                cells = sum(G.faces.neighbors(faces, :), 2);
                comareas = mf(:,3);
                ratio = comareas./G.faces.areas(faces);
                hT = T(faces) .* ratio;
                nnc.cells(:, j) = cells;
                nnc.hT(:, j)    = hT;
            end
            nnc.T  = 1./(1./nnc.hT(:,1) + 1./nnc.hT(:,2));
        end
        
        function intXn = computeIntxnRelation(nwm)
            % Compute intersection relations between subgrids
            fprintf([' -- Computing intersection relations between ',...
                'subgrids\n'])
            % Non-mactching faces
            nmf_CV = [nwm.nonMatchingIntxnRelation([1,2], 'top'); ...
                nwm.nonMatchingIntxnRelation([1,2], 'bot')];
            nmf_CV = nwm.mapIntxnRelationCV(nmf_CV);

            nmf_VW = [nwm.nonMatchingIntxnRelation([2,3], 'heel'); ...
                nwm.nonMatchingIntxnRelation([2,3], 'toe')];
            nmf_VW = nwm.mapIntxnRelationVW(nmf_VW);

            % Mactching faces
            mf_CV = nwm.matchingIntxnRelation([1,2]);
            mf_CV = nwm.mapIntxnRelationCV(mf_CV);

            mf_VW = nwm.matchingIntxnRelation([2,3]);
            mf_VW = nwm.mapIntxnRelationVW(mf_VW);

            intXn.nonMatchingFaces = [nmf_CV; nmf_VW];
            intXn.matchingFaces = [mf_CV; mf_VW];
        end

        function nmf = nonMatchingIntxnRelation(nwm, grdInd, bdyLoc)
            % Compute intersection relations of non-matching faces on
            % the boundaries of input subgrids
            % gridInd: Indices of the subgrids involving in the computation
            %          [i-1, i], i <= number of subgrids
            % bdyLoc:  The location of boundaries
            %          'top'  | 'bot' for CPG and VOI grid
            %          'heel' | 'toe' for VOI grid and HW grid
            %
            subG = nwm.assignInputSubGrds();
            assert(grdInd(2)<=numel(subG));
            assert(grdInd(2)==grdInd(1)+1);
            [G1, G2] = deal( subG{grdInd(1)}, subG{grdInd(2)} );
            types = nwm.assignInputSubGrdTypes();
            F = G2.parentInfo.faces;
            switch bdyLoc
                case {'top', 'heel'}
                    f1 = F{1};
                    f2 = find(G2.faces.surfaces == 1);
                case {'bot',  'toe'}
                    f1 = F{end};
                    f2 = find(G2.faces.surfaces == max(G2.faces.surfaces));
                otherwise
                    error('Unknow boundary location')
            end
            t1 = clock;
            fprintf('      %8s - %8s: %7s boundary, ', types{grdInd(1)},...
                types{grdInd(2)}, bdyLoc)
            nmf = handleNonMatchingFaces(G1, f1, G2, f2, ...
                'isfaceNodesSorted', true);
            t2 = clock;
            fprintf('elapsed time %.2f [s]\n', etime(t2, t1))
        end

        function mf  = matchingIntxnRelation(nwm, grdInd)
            % Compute intersection relations of matching faces on the
            % layered boundaries of input subgrids
            % gridInd: Indices of the subgrids involving in the computation
            %          [i-1, i], i <= number of subgrids
            
            subG = nwm.assignInputSubGrds();
            assert(grdInd(2)<=numel(subG));
            assert(grdInd(2)==grdInd(1)+1);
            [G1, G2] = deal( subG{grdInd(1)}, subG{grdInd(2)} );
            types = nwm.assignInputSubGrdTypes();
            C  = G2.parentInfo.cells;
            BN = G2.parentInfo.bdyNodes;
            t1 = clock;
            fprintf('      %8s - %8s: layered boundary, ', types{grdInd(1)},...
                types{grdInd(2)})
            mf = handleMatchingFaces(G1, C, BN, G2);
            t2 = clock;
            fprintf('elapsed time %.2f [s]\n', etime(t2, t1))
        end

        function model = setupSimModel(nwm, rock, T_all, N_all)
            % Setup simulation model passed to ad-blackoil simulator for 
            % the global grid
            % rock:  Rock of global grid
            % T_all: Full transmissibility
            % N_all: Neighborship of all connections
            gravity reset on
            G    = nwm.gloGrid;
            f    = nwm.fluid;
            % Internal connections
            intCon = all(N_all, 2);
            N = N_all(intCon, :);
            T = T_all(intCon);
            % Phase component
            ph = nwm.getPhaseFromDeck();
            model = GenericBlackOilModel(G, rock, f, 'water', ph.wat, ...
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
        
        function [hasAquifer, output] = handleAquifers(nwm)
            % Handle the aquifers: only support the Fetkovich aquifers
            deck = nwm.inputDeck;
            hasAquifer = all( isfield(deck.SOLUTION, {'AQUANCON', 'AQUFETP'}) );
            if ~hasAquifer
                output = nan;
                return
            end
            [GC, ~, ~] = nwm.assignInputSubGrds ();
            output = processAquifer(deck, GC);
            % Map the conn cells
            aquifers = output.aquifers;
            aquind = output.aquind;
            cells = aquifers(:, aquind.conn);
            mapc = nwm.cellMapFromInputSubGrdsToGloGrd();
            mapc = mapc{1};
            cells = arrayfunUniOut(@(c)mapc(mapc(:,1) == c, 2), cells);
            idx = ~cellfun(@isempty, cells);
            cells = cell2mat( cells(idx) );
            aquifers = aquifers(idx, :);
            aquifers(:, aquind.conn) = cells;
            output.aquifers = aquifers;
            % The aquifer connects to the VOI grid
            if ~all(idx)
                output = nwm.getAquifersVOIG(output);
            end
            
        end
        
        function output = getAquifersVOIG(nwm, output)
            % Aquifers connected to the VOI grid.
            % Note that the implementation only support the influx bottom aquifer
            % Assign data
            G = nwm.gloGrid;
            [~, GV] = nwm.assignInputSubGrds();
            aquifers = output.aquifers;
            aquind = output.aquind;
            % Get the connection faces and cells (GV)
            nSurf = GV.layers.num+1;
            surfInd = G.faces.surfaces(G.faces.grdID==2);
            assert(nSurf == max(surfInd));
            facesV = find( G.faces.surfaces==nSurf & G.faces.grdID==2 );
            N = G.faces.neighbors(facesV,:);
            assert( all( ~all(N,2) ) );
            connV = sum(N, 2);
            % The connection faces and cells (GC)
            connC = aquifers(:, aquind.conn);
            facesC = zeros(size(connC));
            for i = 1 : numel(connC)
                facePos = G.cells.facePos(connC(i)) : ...
                    G.cells.facePos(connC(i)+1)-1;
                f = G.cells.faces(facePos, :);
                % Note only support the influx bottom aquifer
                facesC(i) = f(f(:,2)==6);
            end
            N = G.faces.neighbors(facesC,:);
            assert( all( ~all(N,2) ) );
            % Get the aquifer alhpa
            deck = nwm.inputDeck;
            aquancon = deck.SOLUTION.AQUANCON;
            influxcoef = cell2mat(aquancon(:, 9));
            influxmultcoef = cell2mat(aquancon(:, 10));
            % Use area weighted (aquifer influx coefficient multiplier=1)
            assert(isnan(influxcoef) & influxmultcoef==-1)
            facesA = [facesC; facesV];
            influxcoef = G.faces.areas(facesA);
            alpha = influxcoef./sum(influxcoef);
            % Assemble the aquifer
            aquifersV = nan(numel(connV), 7);
            aquifersV(:, aquind.conn) = connV;
            aquifersV(:, aquind.depthconn) = G.cells.centroids(connV, 3);
            flds = {'aquid', 'pvttbl', 'J', 'C', 'depthaq'};
            for i = 1 : numel(flds)
                fld = flds{i};
                aquifersV(:, aquind.(fld)) = unique(aquifers(:, aquind.(fld)));
            end
            aquifers = [aquifers; aquifersV];
            aquifers(:, aquind.alpha) = alpha;
            output.aquifers = aquifers;
            output.connFaces = [facesC; facesV];
        end
        
        function schedule = getSimSchedule(nwm, model, varargin)
            % Get the simulation schedule for the global grid from the 
            % production/injection control data in deck 
            fprintf(' -- Converting schedule from input deck\n')
            opt = struct('refDepthFrom', 'deck');
            opt = merge_options(opt, varargin{:});
            G = nwm.gloGrid;
            % Assign the CPG schedule first, to get the well structure of
            % CPG model
            modelC = nwm.setupCPGSimModel();
            scheduleC = nwm.getCPGSimSchedule(modelC);
            [wc, ~, WI] = nwm.getWellCellPara(model);
            % Need to map the well cells of other wells from input CPG to
            % global grid
            mapc = nwm.cellMapFromInputSubGrdsToGloGrd();
            mapc = mapc{1};
            % Define a tmp rock
            rockTmp.perm = nan(G.cells.num, 3);
            for i = 1 : numel(scheduleC.control)
                W0 = scheduleC.control(i).W;
                ii = arrayfun(@(w)strcmp(nwm.well.name, w.name), W0);
                % Map the cells of the other wells
                WRegular = W0(~ii);
                for w = 1 : numel(WRegular)
                    c0 = WRegular(w).cells;
                    [tmp, ~, ic] = intersect(c0, mapc(:,1), 'stable');
                    assert(all(tmp==c0));
                    WRegular(w).cells = mapc(ic, 2);
                end
                % Well structure for the HW
                switch opt.refDepthFrom
                    case 'deck'
                        refDepth = W0(ii).refDepth;
                        fprintf(['    Info : The reference depth of %s',...
                            ' adopts the value from deck\n'], nwm.well.name)
                    case 'trajectory'
                        refDepth = nwm.well.trajectory(1,3);
                        fprintf(['    Info : The reference depth of %s',...
                            ' has been set to the depth of first well',...
                            ' point\n'], nwm.well.name)
                    case 'topNode'
                        refDepth = nwm.wellboreGrid.cells.centroids(1,3);
                        fprintf(['    Info : The reference depth of %s',...
                            ' has been set to the depth of top node',...
                            ' in multi-segment well definition\n'], ...
                            nwm.well.name)
                    otherwise
                        error('Unknown reference depth definition type')
                end
                % Redefine some fields
                W             = W0(ii);
                W.cells       = wc;
                [W.r, W.rR]   = deal( nan(numel(wc),1) );
                W.dir         = repmat('X', numel(wc), 1);
                W.WI          = WI;
                W.refDepth    = refDepth;
                W.cell_origin = ones(numel(wc),1);
                W.cstatus     = true(numel(wc),1);
                % Call 'addWell' to compute the 'W.dz'
                WTmp = addWell([], G, rockTmp, wc, 'name', 'WTmp', ...
                    'refDepth', refDepth);
                W.dZ = WTmp.dZ;
                % Combine with the other wells
                WNew = [W; WRegular];
                scheduleC.control(i).W = WNew;
            end
            schedule = scheduleC;
        end
        
        function model = setupCPGSimModel(nwm)
            % Setup simulation model passed to ad-blackoil simulator for
            % the input CPG grid
            f    = nwm.fluid;
            [GC, ~] = nwm.assignInputSubGrds();
            rockC = nwm.getCPGRockFromDeck();
            % Phase component
            ph = nwm.getPhaseFromDeck();
            model = GenericBlackOilModel(GC, rockC, f, 'water', ph.wat, ...
                'oil', ph.oil, 'gas', ph.gas, 'vapoil', ph.vapo, 'disgas', ph.disg);
        end
        
        function schedule = getCPGSimSchedule(nwm, model)
            % Get the simulation schedule for input CPG from deck
            deck = nwm.inputDeck;
            schedule = convertDeckScheduleToMRST(model, deck);
        end
        
        function [wc, wf, WI] = getWellCellPara(nwm, model)
            %  Get parameters for well cells of the HW
            %  wc: well cell indices
            %  wf: well face indices
            %  WI: well indices of well cells
            wc = nwm.getWellCells();
            G = nwm.gloGrid;
            % Wellbore faces (always the first face of wc)
            wf = arrayfun(@(c)G.cells.faces(G.cells.facePos(c),1), wc);
            assert(all( ~all( G.faces.neighbors(wf, :), 2) ))
            % Well index
            WI = model.operators.T_all(wf);
        end
        
        function wc = getWellCells(nwm)
            % Get well cell indices
            G = nwm.gloGrid;
            [~, ~, GW] = nwm.assignSubGrds();
            nA = GW.radDims(1);
            % All HW grid cells at global grid
            mapc = nwm.cellMapFromSubGrdsToGloGrd();
            cells = mapc{3}(:,2);
            assert(all(G.cells.grdID(cells)==3));
            cells = reshape(cells, [], GW.layers.num);
            % Completed segments
            segs = nwm.well.openedSegs;
            assert(all(segs <= nwm.well.segmentNum))
            % Well cells 
            wc = cells(1:nA, segs);
            wc = wc(:);
        end
        
        function state0 = getInitState(nwm, model)
            % Get initial state by equilibrium initialization
            deck = nwm.inputDeck;
            state0 = initStateDeck(model, deck);
        end
        
        function ph = getPhaseFromDeck(nwm)
            % Get phase components from input deck
            deck = nwm.inputDeck;
            rspec = deck.RUNSPEC;
            ph = struct('wat', false, 'oil', false, 'gas', false, ...
                'vapo', false, 'disg', false);
            ph.wat = isfield(rspec, 'WATER') && rspec.WATER;
            ph.oil = isfield(rspec, 'OIL') && rspec.OIL;
            ph.gas = isfield(rspec, 'GAS') && rspec.GAS;
            ph.vapo = isfield(rspec, 'VAPOIL') && rspec.VAPOIL;
            ph.disg = isfield(rspec, 'DISGAS') && rspec.DISGAS;
        end
        
        function mapc = cellMapFromInputSubGrdsToGloGrd(nwm)
            % Cell map from input subgrids to global grid
            subG = nwm.assignSubGrds();
            mapc = nwm.cellMapFromSubGrdsToGloGrd();
            for i = 1 : numel(subG)
                mapc{i}(:,1) = subG{i}.cells.map(mapc{i}(:,1));
            end
        end

        function mapf = faceMapFromInputSubGrdsToGloGrd(nwm)
            % Face map from input subgrids to global grid
            subG = nwm.assignSubGrds();
            mapf = nwm.faceMapFromSubGrdsToGloGrd();
            for i = 1 : numel(subG)
                mapf{i}(:,1) = subG{i}.faces.map(mapf{i}(:,1));
            end
        end

        function mapc = cellMapFromSubGrdsToGloGrd(nwm)
            % Cell map from updated subgrids (after-remove cells) to global
            % grid
            subG = nwm.assignSubGrds();
            nc = cellfun(@(G)G.cells.num, subG);
            if size(nc,1) > size(nc,2); nc = nc'; end
            csnc = cumsum([0, nc]);
            mapc = arrayfunUniOut(@(i)[(1:nc(i))', (1:nc(i))'+csnc(i)], ...
                (1:numel(nc))');
        end

        function mapf = faceMapFromSubGrdsToGloGrd(nwm)
            % Face map from updated subgrids (after-remove cells) to global
            % grid
            subG = nwm.assignSubGrds();
            nf = cellfun(@(G)G.faces.num, subG);
            if size(nf,1) > size(nf,2); nf = nf'; end
            csnf = cumsum([0, nf]);
            mapf = arrayfunUniOut(@(i)[(1:nf(i))', (1:nf(i))'+csnf(i)], ...
                (1:numel(nf))');
        end

        function grdecl = getGrdEclFromDeck(nwm)
            % Get ECLIPSE grid structure from deck
            deck = nwm.inputDeck;
            % Does not assign the 'ACTNUM' for robustness
            if isfield(deck.GRID, 'COORD') % CPG
                fn = {'cartDims', 'COORD', 'ZCORN'};
            else % Cartesian grid
                fn = {'cartDims', 'DX', 'DY', 'DZ', 'TOPS'};
            end
            for i = 1 : numel(fn)
                grdecl.(fn{i}) = deck.GRID.(fn{i});
            end
        end
        
        function checkCellMaps(nwm)
            % Check the cell maps
            G  = nwm.gloGrid;
            % From updated subgrids (after-remove cells) to global grid
            mapc = nwm.cellMapFromSubGrdsToGloGrd();
            subG = nwm.assignSubGrds();
            for i = 1 : numel(subG)
                p1 = subG{i}.cells.centroids(mapc{i}(:,1),:);
                p2 = G.cells.centroids(mapc{i}(:,2),:);
                assert( all( all(p1==p2) ), 'Wrong cell map!')
            end
            % From input subgrids to global grid
            mapc = nwm.cellMapFromInputSubGrdsToGloGrd();
            subG = nwm.assignInputSubGrds();
            for i = 1 : numel(subG)
                p1 = subG{i}.cells.centroids(mapc{i}(:,1),:);
                p2 = G.cells.centroids(mapc{i}(:,2),:);
                assert( all( all(p1==p2) ), 'Wrong cell map!')
            end
            disp('   Corrected cell map   ')
        end

        function checkFaceMaps(nwm)
            % Check the face maps
            G  = nwm.gloGrid;
            % From updated subgrids (after-remove cells) to global grid
            mapf = nwm.faceMapFromSubGrdsToGloGrd();
            subG = nwm.assignSubGrds();
            for i = 1 : numel(subG)
                p1 = subG{i}.faces.centroids(mapf{i}(:,1),:);
                p2 = G.faces.centroids(mapf{i}(:,2),:);
                assert( all( all(p1==p2) ), 'Wrong face map!')
            end
            % From input subgrids to global grid
            mapf = nwm.faceMapFromInputSubGrdsToGloGrd();
            subG = nwm.assignInputSubGrds();
            for i = 1 : numel(subG)
                p1 = subG{i}.faces.centroids(mapf{i}(:,1),:);
                p2 = G.faces.centroids(mapf{i}(:,2),:);
                assert( all( all(p1==p2) ), 'Wrong face map!')
            end
            disp('   Corrected face map   ')
        end

        function [errf1, errf2] = checkIntxnRelation(nwm, intXn)
            % Check the intersection relation by comparing the face areas
            G  = nwm.gloGrid;
            R  = [intXn.nonMatchingFaces(:, 1:3);intXn.matchingFaces];
            f1 = unique( R(:,1) );
            f2 = unique( R(:,2) );
            A  = R(:,3);
            Af1  = G.faces.areas(f1);
            Af2  = G.faces.areas(f2);
            Af1_ = arrayfun(@(f)sum(A(R(:,1) == f)), f1);
            Af2_ = arrayfun(@(f)sum(A(R(:,2) == f)), f2);
            errf1 = abs(Af1 - Af1_) ./ Af1;
            errf2 = abs(Af2 - Af2_) ./ Af2;

            assert( all( ~all(G.faces.neighbors(R(:,1),:),2) ) )
            assert( all( ~all(G.faces.neighbors(R(:,2),:),2) ) )
        end

        function plotNonMatchingIntxnRelation(nwm, intXn, f1)
            % Plot the intersection relations of non-matching face f1
            G = nwm.gloGrid;
            nmf = intXn.nonMatchingFaces;
            c1  = sum(G.faces.neighbors(f1, :),2);
            idx = nmf(:,1)== f1;
            if nnz(idx) > 0
                cla, hold on
                f2  = nmf(idx, 2);
                c2  = sum(G.faces.neighbors(f2, :),2);
                plotGrid(G, c1, 'facecolor', 'none')
                arrayfun(@(c)plotGrid(G, c, 'facecolor', rand(3,1)), c2)
                cc = nmf(idx, 4:6);
                plot3(cc(:,1), cc(:,2), cc(:,3), 'rs', ...
                    'markersize', 6, 'markerFacecolor', 'k')
                view(3),axis off
            else
                warning(['Face %d does not appear in the non-matching ', ...
                    'intersection relations'], f1)
            end
        end

        function plotMatchingIntxnRelation(nwm, intXn, f1)
            % Plot the intersection relations of matching face f1
            G = nwm.gloGrid;
            mf = intXn.matchingFaces;
            c1  = sum(G.faces.neighbors(f1, :),2);
            idx = mf(:,1)== f1;
            if nnz(idx) > 0
                f2  = mf(idx, 2);
                c2  = sum(G.faces.neighbors(f2, :),2);
                cla, hold on
                plotGrid(G, c1, 'facecolor', 'none')
                arrayfun(@(c)plotGrid(G, c, 'facecolor', rand(3,1)), c2)
                view(3),axis off
            else
                warning(['Face %d does not appear in the matching ', ...
                    'intersection relations'], f1)
            end
        end

    end

    methods (Access = protected)

        function CV = mapIntxnRelationCV(nwm, CV)
            % Map the intersection relation from input subgrids (GC and GV)
            % to global grid
            % CV(:,1): faces of GC
            % CV(:,2): faces of GV
            mapf = nwm.faceMapFromInputSubGrdsToGloGrd();
            [mapfC, mapfV] = deal(mapf{1}, mapf{2});
            fC  = arrayfunUniOut(@(f)mapfC(f == mapfC(:,1), 2), CV(:,1));
            fV1 = arrayfun(@(f)mapfV(f == mapfV(:,1), 2), CV(:,2));
            % When the VOI surfaces constitute part of the top or bottom
            % boundary of the global grid, fC will be empty
            idx = ~cellfun(@isempty, fC);
            fC  = fC(idx);
            if ~isempty(fC)
                fC   = cell2mat(fC);
                fV1  = fV1(idx);
                ia   = CV(idx, 3:end);
                CV   = [fC, fV1, ia];
            else
                CV   = [];
            end
        end

        function VW = mapIntxnRelationVW(nwm, VW)
            % Map the intersection relation from input subgrids (GV and GW)
            % to global grid.
            % VW(:,1): faces of GV
            % VW(:,2): faces of GW
            mapf = nwm.faceMapFromInputSubGrdsToGloGrd();
            [mapfV, mapfW] = deal(mapf{2}, mapf{3});
            fV2  = arrayfun(@(f)mapfV(f == mapfV(:,1), 2),  VW(:,1));
            fW   = arrayfun(@(f)mapfW(f == mapfW(:,1), 2),  VW(:,2));
            ia   = VW(:,3:end);
            VW   = [fV2, fW, ia];
        end

        function GCu = updateCPG(nwm)
            % Get the updated CPG whose cells inside the VOI are removed
            [GC, GV] = nwm.assignInputSubGrds();
            cV = cell2mat(GV.parentInfo.cells);
            [GCu, mapc, mapf] = removeCells(GC, cV);
            % Assign some field to keep consistency of the subgrids
            ijk = gridLogicalIndices(GCu);
            GCu.cells.layers = ijk{3};
            GCu.faces.surfaces = nan(GCu.faces.num,1);
            GCu.cells.map  = mapc;
            GCu.faces.map  = mapf;
        end

        function GVu = updateVOIGrid(nwm)
            % Get the updated VOI grid whose cells inside the HW region are
            % removed
            [~, GV, GW] = nwm.assignInputSubGrds();
            cW = cell2mat(GW.parentInfo.cells);
            [GVu, mapc, mapf] = removeCells(GV, cW);
            % Map the ID of layers and surfaces for cells and faces
            GVu.cells.layers   = GVu.cells.layers(mapc);
            GVu.faces.surfaces = GVu.faces.surfaces(mapf);
            GVu.cells.map  = mapc;
            GVu.faces.map  = mapf;
        end

        function GWu = updateHWGrid(nwm)
            % Get the updated HW grid
            % No cells are removed
            [~, ~, GW] = nwm.assignInputSubGrds();
            GWu = GW;
            % Assign some field to keep consistency of the subgrids
            GWu.cells.map = (1:GW.cells.num)';
            GWu.faces.map = (1:GW.faces.num)';
        end

    end
end

function checkDeck(deck)
    if ~isempty(fieldnames(deck.REGIONS))
        error(['The near-wellbore model now only',...
            ' supports single region definition'])
    end
    if ~isfield(deck.SOLUTION, 'EQUIL')
        error(['The keyword ''EQUIL'' in the deck input is required ',...
            'for the equilibrium initialization'])
    end
    if size(deck.SOLUTION.EQUIL, 1) > 1
        error(['The near-wellbore model now only',...
            ' supports single equilibration region'])
    end
end
