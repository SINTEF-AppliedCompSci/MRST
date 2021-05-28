function schedule = upscaleSchedule(model, schedule, varargin)
%Upscale a schedule to a coarser model
%
% SYNOPSIS:
%   schedule = upscaleSchedule(model, schedule)
%
% REQUIRED PARAMETERS:
%   model      - The coarse model the schedule is to be converted to.
%                Assumed to be derived from the fine model used with
%                schedule.
%
%   schedule   - Schedule to be upscaled.
%
% OPTIONAL PARAMETERS:
%   'wellUpscaleMethod' - Upscaling method for well indices. The default is
%                         to recompute the well indices in the new block.
%                         Other options are 'sum', 'harmonic' and 'mean'.
%                         We recommend applying a dedicated upscaling
%                         routine and replacing these values if well-bore
%                         flow performance is important.
%
%   'bcUpscaleMethod'   - Interpolation strategy used for boundary
%                         conditions. Possible options:
%                           - linear: Default.
%                           - idw:    Inverse distance weighting
%                           - mean:   Mean value.
%                           - nearest: Nearest neighbor.
%                         In addition, any unknown arguments will be passed
%                         onto the interpolation routine used. Depending on
%                         the dimensionality and the boundary conditions,
%                         this is either `interp1` or `scatteredInterp`.
%
% RETURNS:
%   schedule   - Schedule upscaled for the coarse model.
%
% NOTE:
%   Support for boundary conditions relies on interpolation. Results should
%   be examined before use for complex grids.
%
% SEE ALSO:
%   `upscaleModelTPFA`, `upscaleState`, `generateCoarseGrid`


%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    opt = struct('wellUpscaleMethod', 'recompute', ...
                 'bcUpscaleMethod',   'linear');
    opt = merge_options(opt, varargin{:});

    for i = 1:numel(schedule.control)
        % Treat wells
        if isfield(schedule.control(i), 'W')
            W = schedule.control(i).W;
            W_coarse = [];
            for j = 1:numel(W)
                w = handleWell(model, W(j), opt);
                W_coarse = [W_coarse; w]; %#ok
            end
            schedule.control(i).W = W_coarse;
        end
        % Treat boundary conditions
        if isfield(schedule.control(i), 'bc')
            bc_coarse = handleBC(model, schedule.control(i).bc, opt);
            schedule.control(i).bc = bc_coarse;
        end
        % Treat source terms
        if isfield(schedule.control(i), 'src')
            src_coarse = handleSRC(model, schedule.control(i).src, opt);
            schedule.control(i).src = src_coarse;
        end
    end
end

function Wc = handleWell(model, W, opt)
    % Handle a single well
    p = model.G.partition;
    Wc = W;
    
    pc = p(W.cells);
    
    % cells
    [Wc.cells, firstInd, newMap] = uniqueStable(pc);
    nc = numel(Wc.cells);
    
    counts = accumarray(newMap, 1);
    
    % Take first direction uncritically (it shouldn't really matter when
    % the well index has been upscaled).
    Wc.dir = W.dir(firstInd);
    if numel(W.r) > 1
        Wc.r = W.r(firstInd);
    end
    % Upscale well index
    if strcmpi(opt.wellUpscaleMethod, 'recompute')
        Wc.WI = computeWellIndex(model.G, model.rock, Wc.r, Wc.cells, 'Dir', Wc.dir);
    else
        switch lower(opt.wellUpscaleMethod)
            case 'sum'
                fn = @(WI, map, counts) accumarray(map, WI);
            case 'harmonic'
                fn = @(WI, map, counts) 1./(accumarray(map, 1./WI)./counts);
            case 'mean'
                fn = @(WI, map, counts)accumarray(map, WI)./counts;
                otherwise
            error(['Unknown upscale mode: "', opt.wellUpscaleMethod, '"'])
        end
        Wc.WI = fn(W.WI, newMap, counts);
    end
    % dZ
    if model.G.griddim > 2
        z = model.G.cells.centroids(Wc.cells, 3);
    else
        z = zeros(numel(Wc.cells), 1);
    end
    Wc.dZ = z - W.refDepth;
    
    % cstatus
    Wc.cstatus = true(nc, 1);
    % Mapping -> coarse to fine
    Wc.fperf = newMap;
    
    if isfield(W, 'topo')
        % Extract topology.
        mp = [0; newMap];
        newtopo = mp(W.topo + 1);    
        % eliminate redundant connections due to cell collapsing in coarser
        % model
        newtopo = sort(newtopo, 2);
        newtopo = uniqueStable(newtopo, 'rows');
        newtopo = newtopo(newtopo(:,1) ~= newtopo(:, 2), :);

        newtopo = sortrows(newtopo);

        Wc.topo = newtopo;
    end
    Wc.parentIndices = firstInd;
end

function bc_coarse = handleBC(model, bc, opt)
    bc_coarse = [];
    if isempty(bc)
        return
    end
    CG = model.G;
    G = CG.parent;    
    % Coarse face -> fine face map
    connCoarse = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2) .';
    
    isFaceBC = false(G.faces.num, 1);
    isFaceBC(bc.face) = true;
    
    coarseFaceNo = zeros(G.faces.num, 1);
    coarseFaceNo(CG.faces.fconn) = connCoarse;
    
    coarseFacesBC = unique(connCoarse(isFaceBC(CG.faces.fconn)));
    for i = 1:numel(coarseFacesBC)
        cf = coarseFacesBC(i);
        
        isCurrentCoarse = false(G.faces.num, 1);
        isCurrentCoarse(coarseFaceNo == cf) = true;
        
        % subset of boundary condition corresponding to current face
        act = isCurrentCoarse(bc.face);
        faces = bc.face(act);
        areas = G.faces.areas(faces);
        
        types = bc.type(act);
        type = types{1};
        values = bc.value(act);
        
        sat = [];
        comp = [];
        assert(all(strcmpi(type, types)), ...
            ['Mixture of boundary condition types on the same coarse face.',...
            ' Not possible to upscale.']);
        
        switch lower(type)
            case 'flux'
                if ~isempty(bc.sat)
                    % Flux-weight saturations
                    sat = bsxfun(@times, bc.sat(act, :), values);
                    sat = sum(sat, 1);
                    sat = sat./sum(sat, 2);
                end
                if isfield(bc, 'components')
                    comp = bsxfun(@times, bc.components(act, :), values);
                    comp = sum(comp, 1);
                    comp = comp./sum(comp, 2);
                end
                
                val = sum(values);
            case 'pressure'
                if ~isempty(bc.sat)
                    % Area weight saturations
                    sat = bsxfun(@times, bc.sat(act, :), areas);
                    sat = sum(sat, 1);
                    sat = sat./sum(sat, 2);
                end
                
                if isfield(bc, 'components')
                    comp = bsxfun(@times, bc.components(act, :), areas);
                    comp = sum(comp, 1);
                    comp = comp./sum(comp, 2);
                end
                
                if isfield(G, 'nodes')
                    C = G.nodes.coords;
                else
                    C = G.faces.centroids;
                end
                
                rescale = @(x) bsxfun(@rdivide, bsxfun(@minus, x, min(C)), max(C) - min(C));
                xq = rescale(CG.faces.centroids(cf, :));
                x = rescale(G.faces.centroids(faces, :));
                if numel(faces) == 1
                    % Single value. Coarse is equal to fine.
                    val = values;
                else
                    switch lower(opt.bcUpscaleMethod)
                        case 'idw'
                            % Inverse distance weighted pressure
                            val = interpolateIDW(x, values, xq, 2);
                        case 'mean'
                            % area weighted mean value
                            val = mean(values.*areas)./sum(areas);
                        case 'nearest'
                            % Nearest neighbor interpolation
                            val = nearestNeighbor(x, values, xq);
                        otherwise
                            % We just assume that this is a valid method for
                            % matlabs unstructured interpolation, after
                            % reducing to the actual dimension
                            keep = true(1, G.griddim);
                            for j = 1:G.griddim
                                % Remove degenerate dimensions
                                if all(abs(x(:, j) - xq(j)) <= 2*max(eps(x(:, j)), eps(xq(j))))
                                    keep(j) = false;
                                end
                            end
                            if nnz(keep) > 1
                                I = scatteredInterpolant(x(:, keep), values,...
                                    opt.bcUpscaleMethod, opt.bcUpscaleMethod);
                                val = I(xq(keep));
                            else
                                val = interp1(x(:, keep), values, xq(keep), ...
                                    opt.bcUpscaleMethod, 'extrap');
                            end
                    end
                    if isempty(val) || any(~isfinite(val))
                        % Always works
                        val = nearestNeighbor(x, values, xq);
                    end
                end
            otherwise
                error(['Unable to upscale boundary condition type "', type, '"']);
        end
        
        bc_coarse = addBC(bc_coarse, cf, type, val, 'sat', sat);
        if isfield(bc, 'components')
            if ~isfield(bc_coarse, 'components')
                bc_coarse.components = [];
            end
            bc_coarse.components = [bc_coarse.components; comp];
        end
    end
end

function src_coarse = handleSRC(model, src, opt)
    src_coarse = [];
    if isempty(src)
        return
    end
    CG = model.G;

    cells = src.cell;
    cpart = CG.partition(cells);
    coarsecells = unique(cpart);
    for i = 1:numel(coarsecells)
        cc = coarsecells(i);
        act = cpart == cc;
        values = src.rate(act);
        assert(all(values < 0) || all(values >= 0), ...
            ['Producing and injecting source terms in same coarse block.',...
            ' Unable to upscale.']);
        val = sum(values);
        if isempty(src.sat)
            sat = [];
        else
            sat = sum(bsxfun(@times, src.sat(act, :), values), 1);
            sat = sat./sum(values);
        end
        % Add coarse source
        src_coarse = addSource(src_coarse, cc, val, 'sat', sat);
    end
end


function val = nearestNeighbor(x, value, xq)
    dist = sqrt(sum(bsxfun(@minus, x, xq).^2, 2));
    [~, minIndex] = min(dist);
    val = value(minIndex);
end