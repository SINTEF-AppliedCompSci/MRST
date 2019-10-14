function limiter = getLimiter(model, type, varargin)
    
    opt = struct('tol', 0, 'limits', []);
    opt = merge_options(opt, varargin{:});
    switch type
        case 'tvb'
            interpSetup = getInterpolationSetup(model, opt);
            limiter = @(state, name, tol, limits) tvb(state, model, interpSetup, name, tol, limits);
        case 'scale'
            limiter = @(state, name, tol, limits) scale(state, model, name, tol, limits);
    end
    
end

%-------------------------------------------------------------------------%
% Common helpers
%-------------------------------------------------------------------------%
function state = update(model, state, name, dof)

    state = model.setProp(state, [name, 'dof'], dof);
    m     = model.discretization.getCellMean(state, dof);
    state = model.setProp(state, name, m);
    
end

%-------------------------------------------------------------------------%
function plotLimiter(disc, state, state0, type)
    if ishandle(1)
        set(0, 'CurrentFigure', 1);
    else
        figure(1);
    end
    clf; hold on
    plotSaturationDG(disc, state0, 'n', 500, 'plot1d', true, 'color', 'k', 'linew', 2);
    plotSaturationDG(disc, state, 'n', 500, 'plot1d', true, 'color', 'k', 'linew', 4, 'LineStyle', '--'); hold off
    ylim([-0.2, 1.2]);
    legend({['Before ', type], ['After ', type]});
    box on;
    drawnow
end
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% TVB limiter
%-------------------------------------------------------------------------%
function state = tvb(state, model, interpSetup, name, tol, limits)

    dof = model.getProp(state, [name, 'dof']);
    disc = model.discretization;
    dim = min(disc.G.griddim, disc.basis.nDof-1);
    state0 = state;
    
    nc = size(dof,2);
    bad = false;
    for i = 1:nc
        [jumpVal, ~, cells] = getInterfaceJumps(disc, dof(:, i), state);
        j = accumarray(cells(:), repmat(jumpVal,2,1) > tol) > 0;
        jump = false(disc.G.cells.num,1);
        jump(cells(:))          = j(cells(:));
        jump(all(state.degree == 0,2)) = false;
        bad = bad | jump;
    end

    if any(bad)
        for i = 1:nc
            dofbar = approximatGradient(model, interpSetup, state, dof(:,i));
            dofbar = dofbar(bad,:)';
            dofbar = dofbar(:);
            dofbar(isnan(dofbar)) = 0;
            ix = disc.getDofIx(state, (1:dim)+1, bad);
            dof(ix,i) = dofbar;
        end
        if any(disc.degree > 1)
            ix = disc.getDofIx(state, (dim+2):disc.basis.nDof, bad);
            dof(ix,:) = 0;
        end
        state = update(model, state, name, dof);
    end
    
%     if opt.plot
%         plotLimiter(disc, state, state0, 'TVB');
%     end
    
    
end

%-------------------------------------------------------------------------%
function [jumpVal, faces, cells] = getInterfaceJumps(disc, dof, state)
    % Get interface jumps for all internal connections in grid

    G     = disc.G;
    faces = find(disc.internalConn);
    cells = G.faces.neighbors(disc.internalConn,:);
    if isfield(G, 'mappings')
        % We assume we are using reordering, and thus only check
        % interface jumps for interfaces against already solved
        % cells
        order   = G.mappings.cellMap.localOrder;
        isUpstr = any(order <= order(~G.cells.ghost)',2);
        keep    = all(isUpstr(cells),2);
        faces   = faces(keep);
    end

    % Saturation function
    s = @(x, c) disc.evaluateDGVariable(x, c, state, dof);

    % Get reference coordinates
    xF    = G.faces.centroids(faces,:);            
    cells = G.faces.neighbors(faces,:);                
    cL    = cells(:,1);
    cR    = cells(:,2);

    % Find inteface jumps
    jumpVal = abs(s(xF, cL) - s(xF, cR));

end

%-------------------------------------------------------------------------%
function dofbar = approximatGradient(model, interpSetup, state, dof)

    G    = model.parentModel.G;
    disc = model.discretization;
    useMap = isfield(G, 'mappings');
    if useMap
        map = G.mappings.cellMap;
        G   = G.parent;
    else
        map = struct('keep', true(G.cells.num,1));
    end
    q    = nan(G.cells.num,1);
    ix   = disc.getDofIx(state, 1);
    q(map.keep) = dof(ix);
    
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
     
    for dNo = 1:disc.dim
        b = interpSetup.tri_basis{dNo};
        for l = 1:disc.dim+1
            loc_cells = interpSetup.C(:, l);
            ccl = interpSetup.tri_cells(loc_cells);
            ds = b(:, l).*q(ccl);
            sigma{dNo} = sigma{dNo} + ds;
        end
    end

    isLin = false(G.cells.num,1);
    isLin(map.keep) = any(state.degree > 0,2);

    nLin = min(G.griddim, disc.basis.nDof-1);
    dofbar = zeros(G.cells.num, nLin);
    db     = nan(sum(interpSetup.cell_support_count+1),1);
    
    a   = [0;cumsum(interpSetup.cell_support_count+1)]+1;
    six = mcolon(a(1:end-1), a(2:end)-2);
    dix = cumsum(interpSetup.cell_support_count+1);
    for dNo = 1:nLin
        
        ix = disc.getDofIx(state, 1+dNo, Inf, true);
        d = zeros(G.cells.num,1);
        d(isLin) = dof(ix(ix>0));
        
        s = sigma{dNo}(vertcat(interpSetup.cell_support{:}));
        
        db(six) = s;
        db(dix) = d;
        [dbmin, dbmax, ~, ~, neg, pos] = getMinMax(db, interpSetup.cell_support_count+1);
        dofbar(:, dNo) = neg.*dbmax + pos.*dbmin;
    end
    
    if useMap
        dofbar = dofbar(map.keep, :);
    end

end

%-------------------------------------------------------------------------%
function interpSetup = getInterpolationSetup(model, opt)

%     model.G = model.discretization.G;
    m = model;
    while isempty(m.operators)
        m = model.parentModel;
    end
%     model.operators = opt.operators;
    w = WENOUpwindDiscretization(m, m.G.griddim);
    [C, pts, cells, basis, supports, linear_weights, scaling] = w.getTriangulation(m);

    interpSetup = struct();
    
    interpSetup.tri_cells = cells;
    interpSetup.tri_basis = basis;
    interpSetup.tri_points = pts;
    interpSetup.linear_weights = linear_weights;
    interpSetup.cell_support = supports;
    interpSetup.scaling = scaling;
    interpSetup.C = C;
    interpSetup.cell_support_count = cellfun(@numel, interpSetup.cell_support);
    
end

%-------------------------------------------------------------------------%
% Scale limiter
%-------------------------------------------------------------------------%
function state = scale(state, model, name, tol, limits)

    G = model.G;
    disc = model.discretization;

    dof = model.getProp(state, [name, 'dof']);
    v   = model.getProp(state, name);
    
    nc = size(dof,2);
    
    for i = 1:nc
        [vMin, vMax] = disc.getMinMax(state, dof(:,nc));
        theta = [(v(:,nc) - limits(:,1))./(v(:,nc) - vMin), ...
                 (limits(:,2) - v(:,nc))./(vMax - v(:,nc)), ...
                 ones(G.cells.num,1)          ];
        theta(~isfinite(theta) | abs(theta) > 1) = 1;
        %             theta(abs(theta) > 1) = 1;
        theta = min(theta, [], 2);

        for dofNo = 2:disc.basis.nDof
            ix = disc.getDofIx(state, dofNo, (1:G.cells.num)', true);
            dof(ix(ix>0),nc) = dof(ix(ix>0),nc).*theta(ix>0);
        end

        vix = any(state.degree > 0);
        dix = disc.getDofIx(state, 1, vix);
        dof(dix,nc) = (dof(dix,nc) - v(vix,nc)).*theta(vix) + v(vix,nc);
    end
    
    state = update(model, state, name, dof);
    [vMin, vMax] = disc.getMinMax(state, dof);
    if any(any(vMin < 0 - eps | vMax > 1 + eps,2)) && 0
        warning('Some values still outside [0,1]')
    end

%     ind = theta < 1;
%     if disc.degree > 1 && 0
%         ix = disc.getDofIx(state, (G.griddim+2):nDofMax, ind);
%         dof(ix,:) = 0;
%     end

end