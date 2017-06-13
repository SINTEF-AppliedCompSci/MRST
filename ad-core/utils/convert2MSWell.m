function w = convert2MSWell(w, varargin)
% simple utility to 'convert' standard well-structure to multi-segment type
opt = struct('cell2node', [], 'connDZ', [], 'nodeDepth', [], 'topo', [], ...
             'segLength', [], 'segRoughness', [], 'segFlowModel', [], ...
             'segType', [], ...
             'segDiam', [], 'G', [], 'vol', []);

opt = merge_options(opt, varargin{:});
nc = numel(w.cells); 

% nodes
w.cell2node = opt.cell2node;
if isempty(w.cell2node) % assume 1-1 cell-node
    w.cell2node = speye(nc);
end
if size(w.cell2node, 2) == 1 % vector form
    w.cell2node = sparse(w.cell2node, (1:numel(w.cells))', 1);
end

cellCount = max(full(sum(w.cell2node,2)), 1);
[nn, nc] = size(w.cell2node);
% matrix for accumulting cells to nodes
%accumMap = full(sparse((1:nc)', w.cell2node, 1, nc, nn));

nodes = struct('depth', nan(nn, 1), 'vol', ones(nn,1), 'dist', nan(nn, 1));
% node-depth relative to ref-depth
if isempty(opt.nodeDepth) % take avg of dZ for cells
    nodes.depth = w.refDepth + (w.cell2node*w.dZ)./cellCount;
else
    nodes.depth = opt.nodeDepth;
end

if ~isempty(opt.vol)
    nodes.vol = opt.vol;
end

% connection-depth relative to ref-depth
w.connDZ = opt.connDZ;
if isempty(w.connDZ) % take same as dZ
    w.connDZ = zeros(size(w.dZ));
end

% topology
w.topo = opt.topo;
if isempty(w.topo) % set to single branch 
    w.topo = [(1:(nn-1))', (2:nn)'];
end
ns = size(w.topo, 1);
%assert(w.topo(1,1)==0, 'Top-segment must be first entry in topo-list')

% create segment structure
s = struct('length', zeros(ns,1), 'roughness', zeros(ns,1), ...
           'diam' ,  zeros(ns,1), 'flowModel', zeros(ns,1), ...
           'topo', w.topo);

       % segment lengths
s.length = opt.segLength;
if isempty(s.length) % use centroids of well-cells
    assert(~isempty(opt.G));
    cc = opt.G.cells.centroids(w.cells, :);
    cn = bsxfun(@rdivide, w.cell2node*cc, cellCount);
    %ds = cn(w.topo(2:end,2), :) - cn(w.topo(2:end,1), :);
    %s.length = zeros(ns,1);
    %s.length(2:end) = sqrt( sum(ds.^2 ,2) );
    ds = cn(w.topo(:,2), :) - cn(w.topo(:,1), :);
    s.length = sqrt( sum(ds.^2 ,2) );
end

s.diam = opt.segDiam;
if isempty(s.diam)
    rc = w.r;
    if numel(rc) == 1
        rc = rc*ones(nc,1);
    end
    rn = (w.cell2node*rc)./cellCount;
    %s.diam = zeros(ns,1);
    %s.diam(2:end) = rn(w.topo(2:end,2), :) + rn(w.topo(2:end,1), :);
    s.diam = rn(w.topo(:,2), :) + rn(w.topo(:,1), :);
end

s.deltaDistance = zeros(size(s.length));
if ~isempty(opt.segType)
    % Check all connections leaving tube
    outFromTube = find(opt.segType == 3);
    for i = 1:numel(outFromTube)
        segNo = outFromTube(i);
        % Find the two neighbor nodes
        neighborNodes = s.topo(segNo, :);
        d = 0;
        n = 0;
        for j = 1:2
            % Connected to two nodes - check which segments connected to
            % these two nodes are the tubing, and compute average length
            % over these segments to assign the valve a fraction of the
            % tubing.
            neigh = neighborNodes(j);
            localSegments = find(any(s.topo == neigh, 2));
            for k = 1:numel(localSegments)
                lseg = localSegments(k);
                if opt.segType(lseg) == 1
                    d = d + s.length(lseg);
                    n = n + 1;
                end
            end
        end
        assert(n > 0)

        s.deltaDistance(segNo) = d/n;
    end
end
if ~isempty(opt.segRoughness)
    s.roughness = opt.roughness;
end

if ~isempty(opt.segFlowModel)
    s.flowModel = opt.flowModel;
end
w.nodes    = nodes;
w.segments = s;
% remove w.topo, this is now under segments
w = rmfield(w, 'topo');

w.isMS = true;
end


    