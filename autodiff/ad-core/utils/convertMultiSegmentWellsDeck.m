function W = convertMultiSegmentWellsDeck(W, G, deck_control, varargin)
% Convert existing parsed well to multisegment well spec from deck input.
% These wells are generally not supported in MRST's solvers but the
% information may still be useful.

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('compact', true);
    opt = merge_options(opt, varargin{:});
    assert(isfield(deck_control, 'COMPSEGS'));
    assert(isfield(deck_control, 'WELSEGS'));
    
    compsegs = deck_control.COMPSEGS;
    welsegs = deck_control.WELSEGS;
    if isempty(compsegs) || isempty(welsegs)
        return
    end

    msnames = compsegs(:, 1);
    wnames = {W.name};
    W(1).isMS = false;
    W(1).topo = [];
    W(1).nodes = [];
    W(1).segments = [];
    W(1).perf_map = [];
    W(1).cells_to_nodes = [];
    nms = numel(msnames);
    for i = 1:nms
        name = msnames{i};
        well = find(strcmp(wnames, name));
        assert(numel(well) == 1);
        assert(strcmp(compsegs{i, 1}, welsegs{i, 1}));
        W(well) = convert(W(well), G, compsegs{i, 2}, welsegs{i, 2}, opt);
        W(well).isMS = true;
    end
end

function w = convert(w, G, compsegs, welsegs, opt)
    top = welsegs.header{1};
    h4 = welsegs.header{4};
    if ~strcmpi(h4, 'INC')
        warning(['Only ''INC'' format supported: WELSEGS entry 4 was ', h4]);
    end
    if top ~= w.refDepth
        warning(['Well ', w.name, ' has top node depth = ', num2str(top), ...
                 ' different from reference depth = ', num2str(w.refDepth)])
    end    
    top_tube = welsegs.header{2};
    top_volume = welsegs.header{3};
    
    segno = vertcat(welsegs.segments{:, 1});
    assert(all(diff(segno) == 1));
    branch = vertcat(welsegs.segments{:, 3});
    nb = max(branch);
    assert(all(unique(branch) == (1:nb)'));

    branches = cell(nb, 1);
    for i = 1:nb
        active = branch == i;
        volumes = [];
        tube_depths = [];
        roughness = [];
        diameters = [];
        z_depths = [];
        number = [];
        lengths = [];
        ws = welsegs.segments(active, :);
        
        parent = ws{1, 4};
        if parent == 1
            local_tube_start = top_tube;
            local_tube_z = top;
        else
            for b = 1:(i-1)
                ix = branches{b}.index == parent;
                if any(ix)
                    local_tube_start = branches{b}.tube_depths(ix);
                    local_tube_z = branches{b}.z_depths(ix);
                    break
                end
            end
        end
        current_tube_depth = local_tube_start;
        current_tube_z = local_tube_z;
        assert(all(diff(vertcat(ws{2:end, 4})) == 1))
        for j = 1:size(ws, 1)
            start = ws{j, 1};
            stop = ws{j, 2};
            % parent = ws{j, 4};

            delta_tube = ws{j, 5};
            delta_z = ws{j, 6};
            diameter = ws{j, 7};
            rough = ws{j, 8};
            vol = ws{j, 9};
            num = (start:stop)';
            n = length(num);
            
            number = [number; num]; %#ok
            volumes = [volumes; repmat(vol, n, 1)]; %#ok
            d = cumsum([current_tube_depth + delta_tube; repmat(delta_tube, n-1)]);
            tube_depths = [tube_depths; d]; %#ok
            zd = cumsum([current_tube_z + delta_z; repmat(delta_z, n-1)]);
            z_depths = [z_depths; zd]; %#ok
            lengths = [lengths; repmat(delta_tube, n)]; %#ok
            diameters = [diameters; repmat(diameter, n)]; %#ok
            roughness = [roughness; repmat(rough, n)]; %#ok

            current_tube_depth = current_tube_depth + delta_tube*n;
            current_tube_z = current_tube_z + delta_z*n;
        end
        branches{i} = struct('index', number, ...
                             'volumes', volumes, ...
                             'parent', parent, ...
                             'tube_depths', tube_depths, ...
                             'roughness', roughness, ...
                             'diameters', diameters, ...
                             'z_depths', z_depths, ...
                             'lengths', lengths);
    end
    % We have computed the branches, next we can map that onto perforations
    cells = compseg_to_cells(G, compsegs);
    a = cells > 0;
    cells = cells(a);
    compsegs = compsegs(a, :);
    if true
        [cells, ia] = unique(cells, 'stable');
        compsegs = compsegs(ia, :);
    end
    branch_no_c = vertcat(compsegs{:, 4});
    perf_map = [];
    for i = 1:nb
        B = branches{i};
        active = branch_no_c == i;
        compsegs_i = compsegs(active, :);
        cells_i = cells(active);
        n = numel(cells_i);
        start_depth = vertcat(compsegs_i{:, 5});
        if false
            stop_depth = vertcat(compsegs_i{:, 6});
            midpoints = (start_depth + stop_depth)/2;
            midpoints_seg = B.tube_depths + B.lengths/2;
        else
            midpoints = start_depth;
            midpoints_seg = B.tube_depths;
        end
        segno = zeros(n, 1);
        for j = 1:n
            % Take the closest point.
            [~, ix] = min(abs(midpoints(j) - midpoints_seg));
            segno(j) = B.index(ix);
        end
        perf_map = [perf_map; cells_i, segno];
    end
    % We might have changed the perforations (repeats etc).
    np = size(perf_map, 1);
    orig_perf = zeros(np, 1);
    for i = 1:np
        pos = find(w.cells == perf_map(i, 1));
        orig_perf(i) = pos;
    end
    w.WI = w.WI(orig_perf);
    w.dZ = w.dZ(orig_perf);
    w.r = w.r(orig_perf);
    w.perf_map = orig_perf;

    w.cstatus = w.cstatus;
    if opt.compact
        % Make into simple table:
        % Nodes (with volumes, depths)
        % Segments (roughness, diameters)
        
        % Topology (nodes to nodes connected via segments).
        N = [];
        n_nodes = 0;
        for b = 1:nb
            % TODO: Include roughness, length here.
            B = branches{b};
            % Connect to parent
            N = [N; B.parent, B.index(1)];
            % Internal connections
            N = [N; B.index(1:end-1), B.index(2:end)];
            n_nodes = max(max(B.index), n_nodes);
        end
        n_seg = size(N, 1);
        volumes = zeros(n_nodes, 1);
        depths = zeros(n_nodes, 1);
        branch = zeros(n_nodes, 1);
        
        volumes(1) = top_volume;
        depths(1) = top;
        for b = 1:nb
            % TODO: Include roughness, length here.
            B = branches{b};
            ix = B.index;
            volumes(ix) = B.volumes;
            depths(ix) = B.z_depths;
            branch(ix) = b;
        end
        % Topology for nodes in ms well
        w.topo = N;
        % Data for each node
        w.nodes = struct('depth', depths, 'vol', volumes, 'branch_id', branch);
        w.segments = struct('roughness', roughness, 'diameter', diameters, 'length', lengths);
        w.cells = perf_map(:, 1);
    else
        w.branches = branches;
    end
    % Mapping between perforations and well nodes
    w.cells_to_nodes = perf_map;
end

function cells = compseg_to_cells(G, compsegs)
    [I, J, K] = gridLogicalIndices(G);
    cI = vertcat(compsegs{:, 1});
    cJ = vertcat(compsegs{:, 2});
    cK = vertcat(compsegs{:, 3});
    cells = zeros(numel(cI), 1);
    for i = 1:numel(cI)
        pos = find(I == cI(i) & J == cJ(i) & K == cK(i));
        if ~isempty(pos)
            cells(i) = pos;
        end
    end
end
