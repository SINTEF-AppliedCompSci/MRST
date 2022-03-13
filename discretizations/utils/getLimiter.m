function limiter = getLimiter(model, type, varargin)
    % Get limiter for dG simulations. Currently only supports TVB and
    % scaling

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

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
    m     = model.discretization.getCellMean(state, Inf, dof);
    state = model.setProp(state, name, m);
    
end

%-------------------------------------------------------------------------%
% TVB limiter
%-------------------------------------------------------------------------%
function state = tvb(state, model, interpSetup, name, tol, limits)

    dof = model.getProp(state, [name, 'dof']);
    v   = model.getProp(state, name);
    disc = model.discretization;
    nLin = nnz(sum(disc.basis.k,2) == 1);
    dim = min(disc.G.griddim, nLin);
    
    nc = size(dof,2);
    isSat = false;
    if strcmpi(name, 's')
        isSat = true;
        nc = nc-1;
    end
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
            dofbar = approximateGradient(model, interpSetup, state, dof(:,i), v(:,i));
            dofbar = dofbar(bad,:)';
            dofbar = dofbar(:);
            dofbar(isnan(dofbar)) = 0;
            ix = disc.getDofIx(state, (1:dim)+1, bad);
            dof(ix,i) = dofbar;
            ix = disc.getDofIx(state, 1, bad);
            dof(ix,i) = v(bad,i);
        end
        if any(disc.degree > 1)
            ix = disc.getDofIx(state, (dim+2):disc.basis.nDof, bad);
            dof(ix,:) = 0;
        end
        if isSat
            fill = disc.getFillSat(state);
            dof(:,end) = fill - sum(dof(:,1:end-1),2);
        end
        state = update(model, state, name, dof);
    end
    
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
function dofbar = approximateGradient(model, interpSetup, state, dof, v)

    G    = model.parentModel.G;
    disc = model.discretization;
    useMap = isfield(G, 'mappings');
    if useMap
        map = G.mappings.cellMap;
        G   = G.parent;
    else
        map = struct('keep', true(G.cells.num,1));
    end
    q = nan(G.cells.num,1);
    q(map.keep) = v;

    
    sigma = cell(1, disc.dim);
    [sigma{:}] = deal(0);
     
    for dNo = 1:disc.dim
        b = interpSetup.tri_basis{dNo};
        for l = 1:disc.dim+1
            loc_cells = interpSetup.C(:, l);
            glob_cells = interpSetup.tri_cells(loc_cells);
            ds = b(:, l).*q(glob_cells);
            sigma{dNo} = sigma{dNo} + ds;
            
        end
    end
    supCount = interpSetup.cell_support_count;
    c = rldecode((1:model.G.cells.num)', supCount);
    % Transform to physical coordinates
    s = 0;
    for dNo = 1:disc.dim
        sigma{dNo} = sigma{dNo}.*interpSetup.scaling.scale(c, dNo);
        s = s + sigma{dNo}.*interpSetup.scaling.mapping(c,:,dNo);
    end
    % Transform to dG reference coordinates
    s = s./(G.cells.dx(c,:)/2);
    sigma = mat2cell(s, size(s,1), ones(size(s,2),1));
    isLin = false(G.cells.num,1);
    isLin(map.keep) = any(state.degree > 0,2);

    nLin = nnz(sum(disc.basis.k,2) == 1);
    k = disc.basis.k;
    dofbar = zeros(G.cells.num, nLin);
    db     = nan(sum(interpSetup.cell_support_count+1),1);
    a   = [0;cumsum(interpSetup.cell_support_count+1)]+1;
    six = mcolon(a(1:end-1), a(2:end)-2);
    dix = cumsum(interpSetup.cell_support_count+1);
    for dNo = 1:G.griddim
        if any(k(:,dNo) == 1 & sum(k,2) == 1)
            dofNo = find(k(:,dNo) == 1 & sum(k,2) == 1);
            ix = disc.getDofIx(state, dofNo, Inf, true);
            d = zeros(G.cells.num,1);
            d(isLin) = dof(ix(ix>0));

            s = sigma{dNo}(vertcat(interpSetup.cell_support{:}));
            db(six) = s;
            db(dix) = d;
            if 1
                % Assume approximated gradients == 0 is due to boundary
                % effects
                dd = rldecode(d, interpSetup.cell_support_count+1, 1);
                db(db==0) = dd(db==0);
            end
            [dbmin, dbmax, ~, ~, neg, pos] = getMinMax(db, interpSetup.cell_support_count+1);
            dofbar(:, dofNo-1) = neg.*dbmax + pos.*dbmin;
        end
    end    
    if useMap
        dofbar = dofbar(map.keep, :);
    end

end

%-------------------------------------------------------------------------%
function interpSetup = getInterpolationSetup(model, opt)

    m = model;
    while isempty(m.operators)
        m = model.parentModel;
    end
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
        [vMin, vMax] = disc.getMinMax(state, dof(:,i));
        theta = [(v(:,i) - limits(:,1))./(v(:,i) - vMin), ...
                 (limits(:,2) - v(:,i))./(vMax - v(:,i)), ...
                 ones(G.cells.num,1)          ];
        theta(~isfinite(theta) | abs(theta) > 1) = 1;
        theta = min(theta, [], 2);

        for dofNo = 2:disc.basis.nDof
            ix = disc.getDofIx(state, dofNo, (1:G.cells.num)', true);
            dof(ix(ix>0),i) = dof(ix(ix>0),i).*theta(ix>0);
        end

        vix = state.degree > 0;
        dix = disc.getDofIx(state, 1, vix);
        dof(dix,i) = (dof(dix,i) - v(vix,i)).*theta(vix) + v(vix,i);
    end
    
    state = update(model, state, name, dof);
    [vMin, vMax] = disc.getMinMax(state, dof);
    if any(any(vMin < 0 - eps | vMax > 1 + eps,2)) && 0
        warning('Some values still outside [0,1]')
    end
    
end
