classdef SimpleWell < PhysicalModel
    properties
        W
        allowCrossflow
        allowSignChange
        allowControlSwitching
        
        dpMaxRel
        dpMaxAbs
        dsMaxAbs
    end
    
    methods
        function well = SimpleWell(W, varargin)
            well = well@PhysicalModel([]);
            well.W = W;
            well.allowCrossflow = true;
            well.allowSignChange = false;
            well.allowControlSwitching = true;
            
            well.dpMaxRel = inf;
            well.dpMaxAbs = inf;
            well.dsMaxAbs = inf;
            if nargin > 1
                well = merge_options(well, varargin{:});
            end
        end
        
        function well = updateWell(well, W)
            well.W = W;
        end
        
        function wsol = validateWellSol(well, resmodel, wsol)
            
        end
        
        function counts = getVariableCounts(wm, fld)
            try
                fn = wm.getVariableField(fld);
            catch
                fn = [];
            end
            if isempty(fn)
                counts = 0;
            else
                counts = 1;
            end
        end
        
        function names = getExtraPrimaryVariableNames(well, resmodel)
            names = {};
        end
        
        function [names, types] = getExtraEquationNames(well, resmodel)
            [names, types] = deal({});
        end
        
        function [vars, names] = getExtraPrimaryVariables(well, wellSol, resmodel)
            names = well.getExtraPrimaryVariableNames(resmodel);
            vars = cell(size(names));
            [vars{:}] = well.getProps(wellSol, names{:});
        end
        
        function [weqs, ctrlEq, extra, extraNames, qMass, qVol, wellSol] = computeWellEquations(well, wellSol0, wellSol, resmodel, q_s, bh, varw, pw, mobw, rhow, compw, dt, iteration)
            [weqs, qMass, mix_s, status, cstatus, qVol] = computeWellContributionsSingleWell(well, wellSol, resmodel, q_s, bh, varw, pw, mobw, rhow, compw);
            ctrlEq =  setupWellControlEquationsSingleWell(wellSol, bh, q_s, status, mix_s, resmodel);
            
            % Update well properties which are not primary variables
            toDouble = @(x)cellfun(@double, x, 'UniformOutput', false);
            cq_sDb = cell2mat(toDouble(qMass));
            
            wellSol.cqs     = bsxfun(@rdivide, cq_sDb, resmodel.getSurfaceDensities);
            wellSol.cstatus = cstatus;
            wellSol.status  = status;
            extra = {};
            extraNames = {};
        end
        
        function [names, types] = getWellEquationNames(well, resmodel)
            act = resmodel.getActivePhases();
            names = {'waterWells', 'oilWells', 'gasWells'};
            types = {'perf', 'perf', 'perf'};
            names = names(act);
            types = types(act);
        end
        
        function wellSol = updateConnectionPressureDrop(well, wellSol0, wellSol, model, q_s, bhp, wellvars, p, mob, rho, comp, dt, iteration)
            if iteration ~= 1
                return
            end
            toDb  = @(x)cellfun(@double, x, 'UniformOutput', false);
            rho     = cell2mat(toDb(rho));
            
            
            active = model.getActivePhases();
            numPh = nnz(active);
            rhoS = [model.fluid.rhoWS, model.fluid.rhoOS, model.fluid.rhoGS];
            rhoS = rhoS(active);
            
            b = bsxfun(@rdivide, rho, rhoS);
            w = well.W;
            if ~isfield(w, 'topo')
                nperf = numel(w.cells);
                w.topo = [(0:(nperf-1))', (1:nperf)'];
            end
            qs = wellSol.cqs; %volumetric in-flux at standard conds


            C = wb2in(w);            % mapping wb-flux to in-flux
            wbqs  = abs(C\qs);       % solve to get well-bore fluxes at surface conds
            wbqst = sum(wbqs, 2);   % total wb-flux at std conds
            % if flux is zero - just use compi
            zi = wbqst == 0;
            if any( zi )
                wbqs(zi,:)  = ones(nnz(zi),1)*w.compi;
                wbqst(zi,:) = sum(wbqs(zi,:), 2);
            end
            % Compute mixture at std conds:
            mixs = wbqs ./ (wbqst*ones(1,numPh));
            % compute volume ratio Vr/Vs
            volRat = compVolRat(mixs, p, b, model);
            % Mixture density at connection conds (by using static b's)
            rhoMix = (mixs*rhoS')./volRat;
            % rhoMix is now density between neighboring segments given by
            % topo(:,1)-topo(:,2) computed by using conditions in well-cell
            % topo(:,2). This is probably sufficiently accurate.

            % get dz between segment nodes and bh-node1
            dpt = [0; w.dZ];
            dz  = diff(dpt);
            g   = norm(gravity);
            ddp = g*rhoMix.*dz; % p-diff between connection neighbors
            % well topology assumes we can traverse from top down, but add a loop
            % just in case crazy ordering.
            cdp    = nan(size(ddp));
            cdp(1) = ddp(1);
            its = 0; maxIts = 100;
            while and(any(isnan(cdp)), its<maxIts)
                its = its +1;
                for cnr = 2:numel(cdp)
                    cdp(w.topo(cnr,2)) = cdp(w.topo(cnr,1)) + ddp(cnr);
                end
            end
            if its == maxIts
                error(['Problem with topology for well: ', wellSol.name, '. Segments appear not to be connected'])
            end
            wellSol.cdp = cdp;
        end
        
        function [q_s, bhp, wellSol, withinLimits] = updateLimits(well, wellSol0, wellSol, model, q_s, bhp, wellvars, p, mob, rho, comp, dt, iteration)
            if ~well.allowControlSwitching
                % We cannot change controls, so we return
                return
            end
            if isfield(well.W, 'status') && ~well.W.status
                % Well is inactive
                return
            end
            lims = well.W.lims;
            qs_double = cellfun(@double, q_s);
            qs_t = sum(qs_double);
            
            actPh = model.getActivePhases();
            gasIx = sum(actPh(1:3));
            oilIx = sum(actPh(1:2));
            watIx = 1;
            
            if ~isnumeric(lims)
                lims.vrat = -inf;
                
                if wellSol.sign > 0
                    modes   = {'bhp', 'rate', 'rate'};
                    flags = [double(bhp) > lims.bhp, ...
                             qs_t > lims.rate, ...
                             qs_t < lims.vrat];
                else
                    modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};

                    % insert dummy limits for missing fields
                    missing_fields = {modes{~cellfun(@(x) isfield(lims, x), modes)}};
                    for f = missing_fields
                       lims = setfield(lims, f{:}, -inf);
                    end
                    
                    
                    
                    [q_w, q_o, q_g] = deal(0);
                    if model.water
                        q_w = qs_double(watIx);
                    end
                    
                    if model.oil
                        q_o = qs_double(oilIx);
                    end
                    
                    if model.gas
                        q_g = qs_double(gasIx);
                    end
                    
                    flags = [double(bhp) < lims.bhp,  ...
                        q_o          < lims.orat, ...
                        q_w + q_o    < lims.lrat, ...
                        q_g          < lims.grat, ...
                        q_w          < lims.wrat, ...
                        qs_t         > -lims.vrat];
                end
            else
                modes = {};
                flags = false;
                assert(isinf(lims))
            end
            % limits we need to check (all others than w.type):
            chkInx = ~strcmp(wellSol.type, modes);
            vltInx = find(flags(chkInx), 1);
            if ~isempty(vltInx)
                withinLimits = false;
                modes  = modes(chkInx);
                switchMode = modes{vltInx};
                fprintf('Well %s: Control mode changed from %s to %s.\n', wellSol.name, wellSol.type, switchMode);
                wellSol.type = switchMode;
                wellSol.val  = lims.(switchMode);
            else
                withinLimits = true;
            end
            
            if ~withinLimits
                v  = wellSol.val;
                switch wellSol.type
                    case 'bhp'
                        bhp = assignValue(bhp, v, 1);
                    case 'rate'
                        for ix = 1:numel(q_s)
                            q_s{ix} = assignValue(q_s{ix}, v*well.W.compi(ix), 1);
                        end
                    case 'orat'
                        q_s{oilIx} = assignValue(q_s{oilIx}, v, 1);
                    case 'wrat'
                        q_s{watIx} = assignValue(q_s{watIx}, v, 1);
                    case 'grat'
                        q_s{gasIx} = assignValue(q_s{gasIx}, v, 1);
                end % No good guess for qOs, etc...
            end
        end
        
        function wellSol = updateWellSol(well, wellSol, variables, dx)
            isBHP = strcmpi(variables, 'bhp');
            if any(isBHP)
                dv = dx{isBHP};
                dv = well.limitUpdateRelative(dv, wellSol.bhp, well.dpMaxRel);
                dv = well.limitUpdateAbsolute(dv, well.dpMaxAbs);
                wellSol.bhp = wellSol.bhp + dv;
                variables = variables(~isBHP);
                dx = dx(~isBHP);
            end
            for i = 1:numel(dx)
                wellSol = well.updateStateFromIncrement(wellSol, dx{i}, [], variables{i});
            end
        end
        
        function [wellSol, well_shut] = updateWellSolAfterStep(well, resmodel, wellSol)
            w = well.W;
            % Check if producers are becoming injectors and vice versa. The indexes
            % of such wells are stored in inx.
            wsg = w.sign;
            ssg = sign(getTotalRate(wellSol));
            closed = wsg ~= ssg;
            % A well can be set to zero rate without beeing shut down. We update inx
            % to take into account this fact.
            closed = closed & ~strcmpi(w.type, 'bhp') & w.val ~= 0;
            if closed && ~well.allowSignChange && well.allowControlSwitching
                fprintf('Well %s shut down.\n', w.name);
                wellSol.status = false;
                well_shut = true;
            else
                well_shut = false;
            end
            
            switched = ~strcmpi(wellSol.type, w.type);
            if switched
                fprintf('Well %s has switched from %s to %s.\n', w.name, ...
                                                                 w.type, ...
                                                                 wellSol.type);
            end
        end
        
        function [fn, index] = getVariableField(model, name)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = 1;
            switch(lower(name))
                case 'bhp'
                    fn = 'bhp';
                case 'qos'
                    fn = 'qOs';
                case 'qgs'
                    fn = 'qGs';
                case 'qws'
                    fn = 'qWs';
                otherwise
                    % This will throw an error for us
                    [fn, index] = getVariableField@PhysicalModel(model, name);
            end
        end
        
        function ws = ensureWellSolConsistency(well, ws) %#ok
            % Run after the update step to ensure consistency of variables
        end
    end
end


function C = wb2in(w)
    conn = w.topo(2:end, :);
    % Number of connections between perforations
    nconn = size(conn, 1);
    % Number of perforations
    nperf = numel(w.cells);
    
    if nconn + 1 ~= nperf
        warning(['Mismatch between connection count (', num2str(nconn+1),...
                ') and perforation count (', num2str(nperf), '). Well model', ...
                'Does not appear to be a tree.']);
    end

    id = (1:nperf)';
    % First identity, then honor topology.
    ii = [id; conn(:, 1)];
    jj = [id; conn(:, 2)];
    vv = [ones(nperf, 1); -ones(nconn, 1)]; 

    C = sparse(ii, jj, vv, nperf, nperf);
end

function volRat = compVolRat(mixs, p, b, model)
%
x = mixs;
dg = isprop(model, 'disgas') && model.disgas;
vo = isprop(model, 'vapoil') && model.vapoil;

if dg || vo
    [~, isgas] = model.getVariableField('sg');
    [~, isoil] = model.getVariableField('so');
    
    both = find(isgas | isoil);
    
    g = mixs(:, isgas);
    o = mixs(:, isoil);
    
    if dg
        rsMax = model.fluid.rsSat(double(p));
    else
        rsMax = 0;
    end
    if isa(model, 'ThreePhaseBlackOilModel')
        % Vapoil/disgas
        if vo
            rvMax = model.fluid.rvSat(double(p));
        else
            rvMax = 0;
        end
        
        gor = abs(g./o);
        gor(isnan(gor)) = inf;
        rs = min(rsMax, gor);
        ogr = abs(o./g);
        ogr(isnan(gor)) = inf;
        rv = min(rvMax, ogr);
        d = 1-rs.*rv;
        x(:,isgas) = (x(:,isgas) - rs.*o)./d;
        x(:,isoil) = (x(:,isoil) - rv.*g)./d;
        x(:,both) = x(:,both).*(x(:,both)>0);
    else
        % Only gas dissolution
        x(:,isgas) = x(:,isgas) - rsMax.*o;
        x(:,isgas) = x(:,isgas).*(x(:,isgas)>0);
    end
end

volRat = sum(x./b ,2);
end

function qt = getTotalRate(sol)
   ns = numel(sol);
   qt       = zeros([ns, 1]);
   if ns == 0
       return
   end
   typelist = {'qWs', 'qOs', 'qGs'};
   types    = typelist(isfield(sol(1), typelist));
   for w = 1:ns
      for t = reshape(types, 1, []),
         x = sol(w).(t{1});
         if ~isempty(x),
            qt(w) = qt(w) + x;
         end
      end
   end
end
