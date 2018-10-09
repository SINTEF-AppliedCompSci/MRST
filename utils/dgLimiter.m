function state = dgLimiter(disc, state, bad, type, varargin)

    opt = struct('plot', false);
    opt = merge_options(opt, varargin{:});
    
    nDofMax = disc.basis.nDof;
    G       = disc.G;
    sdof    = state.sdof;
    
    
    if opt.plot
        figure(1); clf; hold on
        plotSaturationDG(disc, state, 'n', 500, 'plot1d', true, 'color', 'b', 'linew', 2);
    end

    switch type
            
         case 'kill'
            % Simple "limiter" that reduces to dG(0) for bad cells
            
            % Fix saturation so that it is inside [0,1]
            ix         = disc.getDofIx(state, 1, bad);
            sdof(ix,:) = min(max(state.s(bad,:), 0), 1);
            % Ensure sum(s,2) = 1
            sdof(ix,:) = sdof(ix,:)./sum(sdof(ix,:),2);
            % Remove dofs for dofNo > 1
            if disc.degree > 0
                ix         = disc.getDofIx(state, 2:nDofMax, bad);
                sdof(ix,:) = [];
                state.degree(bad) = 0;
            end
            % Assign new dofs to state
            state.sdof = sdof;

        case 'tvb'
            % TVB limiter
            
            nPh = size(sdof,2);
            for phNo = 1:nPh
                dofbar = approximatGradient(sdof(:,phNo), state, disc);
                dofbar = dofbar(bad,:)';
                dofbar = dofbar(:);
                dofbar(isnan(dofbar)) = 0;
                ix = disc.getDofIx(state, (1:G.griddim)+1, bad);
                sdof(ix,phNo) = dofbar;
            end
            
            if disc.degree > 1
                ix = disc.getDofIx(state, (G.griddim+2):nDofMax, bad);
                sdof(ix,:) = [];
                state.degree(bad) = 1;
            end
            state.sdof = sdof;
            
        case 'scale'
            
            s = state.s;
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            theta = [(s(:,1) - 0)./(s(:,1) - sMin), ...
                     (1 - s(:,1))./(sMax - s(:,1)), ...
                     ones(G.cells.num,1)          ];
            theta(~isfinite(theta) | abs(theta) > 1) = 1;
%             theta(abs(theta) > 1) = 1;
            theta = min(theta, [], 2);
            
            for dofNo = 2:nDofMax%(1:G.griddim)+1
                ix = disc.getDofIx(state, dofNo, (1:G.cells.num)', true);
                sdof(ix(ix>0),:) = sdof(ix(ix>0),:).*theta(ix>0);
            end
            
            satIx = state.degree > 0;
            dofIx = disc.getDofIx(state, 1, satIx);
            sdof(dofIx,:) = (sdof(dofIx,:) - s(satIx,:)).*theta(satIx) + s(satIx,:);
            state.sdof = sdof;
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            if any(sMin < 0 - eps | sMax > 1 + eps)
                warning('Some saturations still outside [0,1]')
            end
            
            ind = theta < 1;
            if disc.degree > 1 && 0
                ix = disc.getDofIx(state, (G.griddim+2):nDofMax, ind);
                sdof(ix,:) = [];
                state.degree(ind & state.degree > 1) = 1;
%                 sdof(ix,:) = 0;
            end
            state.sdof = sdof;
            state.scaled = ind;
            
        case 'orderReduce'
            
            [sMin, sMax] = disc.getMinMaxSaturation(state);
            outside = sMin < 0 - disc.outTolerance | ...
                      sMax > 1 + disc.outTolerance;
                  
            bad = bad & outside & state.degree > 1;
                  
            if any(bad)
                ix = disc.getDofIx(state, (G.griddim + 2):nDofMax, bad);
                sdof(ix,:) = [];
                state.degree(bad) = 1;
%                 sdof(ix,:) = 0;
                state.sdof = sdof;
            end

    end
    
    % Update state dofPos and satuarion
    state   = disc.updateDofPos(state);
    state.s = disc.getCellSaturation(state);
    
    if opt.plot
        plotSaturationDG(disc, state, 'n', 500, 'plot1d', true, 'color', 'r', 'linew', 4, 'LineStyle', '--'); hold off
        legend({['Before ', type], ['After ', type]});
        drawnow
    end
    
    if 0
        [sum(sdof,2), rldecode((1:G.cells.num)', state.nDof, 1)]
    end
end

function dofbar = approximatGradient(dof, state, disc)

%     ind  = 1:disc.basis.nDof:disc.G.cells.num*disc.basis.nDof

    G    = disc.G;
    useMap = isfield(G, 'mappings');
    if useMap
        map = G.mappings.cellMap;
        G   = G.parent;
    else
        map = struct('keep', true(G.cells.num,1));
    end
    q    = nan(G.cells.num,1);
    dofIx   = disc.getDofIx(state, 1);
    q(map.keep) = dof(dofIx);
    
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
     
    for dNo = 1:disc.dim
        b = disc.interp_setup.tri_basis{dNo};
        for l = 1:disc.dim+1
            loc_cells = disc.interp_setup.C(:, l);
            ccl = disc.interp_setup.tri_cells(loc_cells);
            ds = b(:, l).*q(ccl);
            sigma{dNo} = sigma{dNo} + ds;
        end
    end

    isLin = false(G.cells.num,1);
    isLin(map.keep) = state.degree > 0;
    cellNo = rldecode((1:G.cells.num)', disc.interp_setup.cell_support_count, 1);
    colNo = mcolon(ones(G.cells.num, 1), disc.interp_setup.cell_support_count);
    dofbar = zeros(G.cells.num,G.griddim);
    for dNo = 1:G.griddim
        ix = vertcat(disc.interp_setup.cell_support{:});
        Sigma = sparse(cellNo, colNo, sigma{dNo}(ix));
        dofIx = disc.getDofIx(state, 1+dNo, Inf, true);
        dd = nan(G.cells.num,1);
        dd(isLin) = dof(dofIx(dofIx>0));
        Sigma = full([Sigma, dd]);
        Sigma(Sigma==0) = nan;
        dofbar(:,dNo) = minmod(Sigma);
    end
    
    if useMap
        dofbar = dofbar(map.keep, :);
    end

end

function v = minmod(val)
    valmin = val;
    valmin(isnan(val)) = -1;
    valmax = val;
    valmax(isnan(val)) = 1;
    v = max(val,[],2).*all(valmin < 0,2) + min(val,[],2).*all(valmax > 0,2);
end