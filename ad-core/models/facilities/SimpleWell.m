classdef SimpleWell < PhysicalModel
    properties
        W
        allowCrossflow
        allowSignChange
        allowControlSwitching
    end
    
    methods
        function well = SimpleWell(W, varargin)
            well = well@PhysicalModel([]);
            well.W = W;
            well.allowCrossflow = true;
            well.allowSignChange = false;
            well.allowControlSwitching = true;
            
            if nargin > 1
                well = merge_options(well, varargin{:});
            end
        end
        
        function well = updateWell(well, W)
            well.W = W;
        end
        
        function [vars, names] = getWellPrimaryVariables(well, wellSol, resmodel)
            actPh = resmodel.getActivePhases();
            activeVars = [actPh, true];
            names = {'qWs', 'qOs', 'qGs', 'bhp'};
            names = names(activeVars);
            
            vars = cell(1, nnz(activeVars));
            for i = 1:numel(names)
                if activeVars(i)
                    vars{i} = wellSol.(names{i});
                end
            end
        end
        
        function wellSol = updateConnectionPressureDrop(well, wellSol, model, wellvars, p, sat, rho, comp)
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
            wellSol(k).cdp = cdp;
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
