function [substate0, substate, submodel, subforces, mappings] = getSubsetProblem(model, state, state0, forces, subs, varargin)
    opt = struct('inflowTreatment', 'pressure', ...
                 'outflowTreatment', 'pressure');
    opt = merge_options(opt, varargin{:});
    % Create subproblem defined by subs
    isComp = isa(model, 'ThreePhaseCompositionalModel');
    assert(islogical(subs));
    state = value(state);
    if all(subs)
       substate0 = state0;
       substate = state;
       submodel = model;
       subforces = forces;
       keepWells = (1:numel(forces.W))';
       faceMap = (1:model.G.faces.num)';
       renum = (1:model.G.cells.num)';
    else
        % Coarsen state
        for i = 1:numel(forces.W)
            wc = forces.W(i).cells;
            if any(subs(wc))
                subs(wc) = true;
            end
        end

        [substate, substate0, cellFields] = mapCellFields(model, state, state0, subs);

        % Coarsen model
        op0 = model.operators;
        G = model.G;
        isSeq = isa(model, 'SequentialPressureTransportModel');
        if isSeq
            G = model.pressureModel.G;
            op0 = model.pressureModel.operators;
            if isComp
                eos = model.pressureModel.EOSModel;
            end
            rhoS = model.pressureModel.getSurfaceDensities();
        else
            if isComp
                eos = model.EOSModel;
            end
            rhoS = model.getSurfaceDensities();
        end
        op = op0;

        submodel = model;

        Nc = size(state.pressure, 1);
        keep = false(Nc, 1);
        keep(subs) = true;

        faceList = (1:G.faces.num)';
        ifaceList = faceList(op0.internalConn);

        activeConn = all(keep(op0.N), 2);
        keepFaces = op0.internalConn;
        keepFaces(keepFaces) = activeConn;

        bndKeep = keep(op0.N);
        isBnd = sum(bndKeep, 2) == 1;

        faceMap = [ifaceList(activeConn); ifaceList(isBnd)];

        % Dims
        nf = nnz(activeConn);
        nc = nnz(keep);
        renum = zeros(1, nc);
        renum(keep) = (1:nc);

        T = op0.T(activeConn);
        N = renum(op0.N(activeConn, :));
        M  = sparse((1:nf)'*[1 1], N, .5*ones(nf,2), nf, nc);
        C  = sparse( [(1:nf)'; (1:nf)'], N, ones(nf,1)*[1 -1], nf, nc);

        % Operators
        op.C = C;
        op.M = M;
        upstr =  @(flag, x) faceUpstr(flag, x, N, [nf, nc]);
        op.faceUpstr = upstr;
        if nc == 1
            % Some kind of bug in AD?
            op.Div = @(x) 0;
        else
            op.Div = @(x) C'*x;
        end
        op.Grad = @(x) -C*x;
        op.faceAvg = @(x) M*x;
        op.AccDiv = @(acc, flux) acc + C'*flux;
        op.splitFaceCellValue = @(operators, flag, x) splitFaceCellValue(operators, flag, x, [nf, nc]);
        if isfield(op, 'diag_updated')
            op = rmfield(op, 'diag_updated');
        end

        % Trans
        op.T = T;
        op.T_all = T;
        % Mess with grid
        op.N = N;
        % Mess with other stuff
        op.pv = op0.pv(keep);
        op.internalConn = true(nf, 1);

        G.parent = G;
        G.subset = keep;

        G.cells.centroids = G.cells.centroids(keep, :);
        G.cells.num = nc;

        % Coarsen forces
        subforces = forces;
        V = state.flux(op0.internalConn, :);
        substate.flux = state.flux(keepFaces, :);


        sgn = 1 - 2*keep(op0.N(isBnd, 1));
        volFluxRes = sgn.*V(isBnd, :);

        isProd = sum(volFluxRes, 2) < 0 ;
        isInj = ~isProd;
        
        left_c = op0.N(isBnd, 1);
        right_c = op0.N(isBnd, 2);
        left = ~keep(left_c);
        global_cells = left.*left_c + ~left.*right_c;

        assert(~any(keep(global_cells)))
        % Face average
        propstate = state;
        propstate.s = propstate.s./sum(propstate.s, 2);
        propstate = model.initPropertyContainers(propstate);
        p_phase = model.getProp(propstate, 'PhasePressures');
        p_phase = value(p_phase);
        src_pressure = p_phase(global_cells, :);
        if ~isfield(state, 'rho')
            r = model.getProp(propstate, 'Density');
            state.rho = value(r);
        end
        if ~isfield(state, 'mob')
            m = model.getProp(propstate, 'Mobility');
            state.mob = value(m);
        end
        
        rho = state.rho(global_cells, :);
        sT = sum(state.s(global_cells, :), 2);
        if isfield(state, 'massFlux')
            mf = state.massFlux(model.operators.internalConn, :);
            massFlux = sgn.*mf(isBnd, :);
        else
            massFlux = volFluxRes.*rho.*sT;
        end

        if isComp
            cflux = state.componentFluxes(op0.internalConn, :);
            comp_mass_fraction_flux = cflux(isBnd, :);
            comp_mass_fraction_flux = nan*comp_mass_fraction_flux;
            comp_mass_fraction_flux = comp_mass_fraction_flux./sum(comp_mass_fraction_flux, 2);
            src_comp = eos.getMoleFraction(comp_mass_fraction_flux);
        else
            src_comp = [];
        end
        volFluxSurf = bsxfun(@rdivide, massFlux, rhoS);

        W = forces.W;

        keepWells = false(numel(W), 1);

        renum = ones(Nc, 1);
        renum(keep) = 1:nc;
        for i = 1:numel(W)
            keepw = keep(W(i).cells);
            if any(keepw)
                assert(all(keepw));
                W(i).cells = renum(W(i).cells);
                keepWells(i) = true;
            end
        end
        W = W(keepWells);
        substate.wellSol = substate.wellSol(keepWells);
        substate0.wellSol = substate0.wellSol(keepWells);

        faces = nf + (1:nnz(isBnd))';
        nf = numel(faces);

        src_sat = nan(nf, numel(rhoS));
        bc = struct('face', faces, 'type', [], 'value', [], 'sat', src_sat, 'components', src_comp);

        mob = state.mob(global_cells, :);
        rho = state.rho(global_cells, :);

        bc.mob = mob;
%                 bc.mob = sT.*mob;

        bc.rho = rho;
        bc.pressure = src_pressure;
        nph = size(mob, 2);
        bc.type = cell(1, nf);
        bc.value = zeros(nf, nph);
        bc = setBoundaryConditions(bc, opt.inflowTreatment, isInj, src_pressure, volFluxRes, volFluxSurf);
        bc = setBoundaryConditions(bc, opt.outflowTreatment, isProd, src_pressure, volFluxRes, volFluxSurf);

        if isComp
            bc.x = state.x(global_cells, :);
            bc.y = state.y(global_cells, :);

            bc.xM = sT.*eos.getMassFraction(bc.x);
            bc.yM = sT.*eos.getMassFraction(bc.y);
        end

        isP = strcmpi(bc.type, 'pressure');
        isR = strcmpi(bc.type, 'rflux');

        bc.sat(isP | isR, :) = state.s(global_cells(isP | isR), :);
        % bc.sat = bc.sat./sum(bc.sat, 2);

        bad = any(isnan(bc.sat), 2);
        bc.sat(bad, :) = state.s(global_cells(bad), :);
        if isComp
            bc.components(bad, :) = state.components(global_cells(bad), :);
        end

        op.T_all = vertcat(op.T, op0.T(isBnd));
        op.internalConn = vertcat(op.internalConn, false(nnz(isBnd), 1));

        allfaces = [ifaceList(activeConn); ifaceList(isBnd)];
        substate.flux = state.flux(allfaces, :);
        G.faces.num = numel(op.T_all);
        if sum(isBnd) == 1
            N_bnd = renum(op0.N(isBnd, :))'.*bndKeep(isBnd, :);
        else
            N_bnd = renum(op0.N(isBnd, :)).*bndKeep(isBnd, :);
        end

        G.faces.neighbors = vertcat(op.N, N_bnd);
        G.faces.centroids = [G.faces.centroids(ifaceList(activeConn), :); G.parent.cells.centroids(global_cells, :)];
        G.faces.areas = [G.faces.areas(ifaceList(activeConn)); G.faces.areas(ifaceList(isBnd))];
        subforces.bc = bc;
        subforces.W = W;

        submodel.operators = op;
        submodel.G = G;
        submodel.FluxDiscretization = [];
        submodel.FacilityModel = FacilityModel(submodel);
        submodel.FacilityModel = submodel.FacilityModel.setupWells(subforces.W);
        if isprop(submodel.AutoDiffBackend, 'useMex')
            submodel.AutoDiffBackend.useMex = false;
        end
        submodel = submodel.validateModel();
        
        substate = submodel.reduceState(substate, true);
        substate0 = submodel.reduceState(substate0, true);
    end
    assert(isempty(forces.bc));
    assert(isempty(forces.src));

    mappings = struct('wells', keepWells, ...
                      'cells', subs, ...
                      'boundary', ifaceList(isBnd), ...
                      'cellRenum', renum, ...
                      'faces', faceMap, ...
                      'cellFields', {cellFields});
end

function bc = setBoundaryConditions(bc, option, isType, pressure, volFluxRes, volFluxSurf)
    switch lower(option)
        case 'volume_flux'
            type = 'rflux';
            vf = volFluxRes(isType, :);
            v = sum(vf, 2);
            s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
        case 'mass_flux'
            type = 'flux';
            vf = volFluxSurf(isType, :);
            v = sum(vf, 2);
            s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
        case 'pressure'
            type = 'pressure';
            v = pressure(isType, :);
            vf = volFluxRes(isType, :);
            s = bsxfun(@rdivide, abs(vf), sum(abs(vf), 2));
        case 'component_flux'
            disp('hello_world');
        otherwise
            error('Unknown treatment %s', option);
    end
    [bc.type{isType}] = deal(type);
    if size(v, 2) == 1
        for i = 1:size(v, 2)
            bc.value(isType, i) = v;
        end
    else
        bc.value(isType, :) = v;
    end
    bc.sat(isType, :) = s;
end

function [substate, substate0, cellFields] = mapCellFields(model, state, state0, subs)
    substate = state;
    substate0 = state0;
    
    flds = fieldnames(state);
    is_cell = false(numel(flds), 1);
    for i = 1:numel(flds)
        fn = flds{i};
        d = state.(fn);
        if isnumeric(d) && size(d, 1) == model.G.cells.num
            is_cell(i) = true;
            substate.(fn) = d(subs, :);
            % Also check state0.
            if isfield(substate0, fn)
                substate0.(fn) = state0.(fn)(subs, :);
            end
        end
    end
    cellFields = flds(is_cell);
end