function w = convert2MSWell(w, varargin)
%Utility for Converting Standard Well Structure to Multi-Segment Type
%
% Derives and includes Well Nodes and Well Segments in the well structure.
%
% SYNOPSIS:
%   W = convert2MSWell(W)

%{
#COPYRIGHT#
%}

opt = struct('cell2node'   , [], ...
             'connDZ'      , [], ...
             'nodeDepth'   , [], ...
             'topo'        , [], ...
             'segLength'   , [], ...
             'segRoughness', [], ...
             'segFlowModel', [], ...
             'segType'     , [], ...
             'segDiam'     , [], ...
             'G'           , [], ...
             'vol'         , []);

opt = merge_options(opt, varargin{:});

% nodes
w.cell2node = create_cell2node(numel(w.cells), opt);

cellCount = max(full(sum(w.cell2node, 2)), 1);
[nn, nc]  = size(w.cell2node);

% matrix for accumulating cells to nodes
% accumMap = full(sparse((1:nc)', w.cell2node, 1, nc, nn));

nodes = struct('depth', nan ([nn, 1]), ...
               'vol'  , ones([nn, 1]), ...
               'dist' , nan ([nn, 1]));

% node-depth relative to ref-depth
nodes.depth = node_depth(w, cellCount, opt);

if ~isempty(opt.vol)
    nodes.vol = opt.vol;
end

% connection-depth relative to ref-depth
w.connDZ = opt.connDZ;
if isempty(w.connDZ)
   % Connection depths not defined.  Place at zero depth relative to
   % reference (w.dZ).
   w.connDZ = zeros(size(w.dZ));
end

% create segment structure, including node connectivity.
s = create_segments(node_topology(nn, opt));

% segment geometry
s.length = segment_length  (w, s.topo,     cellCount, opt);
s.diam   = segment_diameter(w, s.topo, nc, cellCount, opt);

s.deltaDistance = segment_deltaDistance(s.length, s.topo, opt);

if ~isempty(opt.segRoughness)
    s.roughness = opt.roughness;
end

if ~isempty(opt.segFlowModel)
    s.flowModel = opt.flowModel;
end

w.nodes    = nodes;
w.segments = s;
w.isMS     = true;
end

%--------------------------------------------------------------------------

function c2n = create_cell2node(nperf, opt)
   c2n = opt.cell2node;

   if isempty(c2n)
      % No pre-defined perf-to-node connectivity.  Assume cell == node.
      c2n = speye(nperf);
   end

   if size(c2n, 2) == 1
      % Perf-to-node given in vector form--expand to matrix form.
      c2n = sparse(c2n, 1 : nperf, 1);
   end
end

%--------------------------------------------------------------------------

function depth = node_depth(w, cellCount, opt)
   if isempty(opt.nodeDepth) % take avg of dZ for cells
      depth = w.refDepth + ((w.cell2node * w.dZ) ./ cellCount);
   else
      depth = opt.nodeDepth;
   end
end

%--------------------------------------------------------------------------

function topo = node_topology(nn, opt)
   topo = opt.topo;

   if isempty(topo) % set to single branch
      topo = [(1 : (nn - 1)).', (2 : nn).'];
   end
end

%--------------------------------------------------------------------------

function s = create_segments(topo)
   ns = size(topo, 1);

   s = struct('length'   , zeros([ns, 1]), ...
              'roughness', zeros([ns, 1]), ...
              'diam'     , zeros([ns, 1]), ...
              'flowModel', zeros([ns, 1]), ...
              'topo'     , topo);
end

%--------------------------------------------------------------------------

function len = segment_length(w, topo, cellCount, opt)
   len = opt.segLength;

   if isempty(len) % use centroids of well-cells
      assert(~isempty(opt.G));

      cc = opt.G.cells.centroids(w.cells, :);
      cn = bsxfun(@rdivide, w.cell2node * cc, cellCount);

      %ds = cn(w.topo(2:end,2), :) - cn(w.topo(2:end,1), :);
      %s.length = zeros(ns,1);
      %s.length(2:end) = sqrt( sum(ds.^2 ,2) );

      ds  = cn(topo(:,2), :) - cn(topo(:,1), :);
      len = sqrt(sum(ds .^ 2, 2));
   end
end

%--------------------------------------------------------------------------

function diam = segment_diameter(w, topo, nc, cellCount, opt)
   diam = opt.segDiam;

   if isempty(diam)
      rc = w.r;

      if numel(rc) == 1
         rc = repmat(rc, [nc, 1]);
      end

      rn = (w.cell2node * rc) ./ cellCount;

      %s.diam = zeros(ns,1);
      %s.diam(2:end) = rn(w.topo(2:end,2), :) + rn(w.topo(2:end,1), :);

      diam = rn(topo(:,2), :) + rn(topo(:,1), :);
   end
end

%--------------------------------------------------------------------------

function dist = segment_deltaDistance(len, topo, opt)
   dist = zeros(size(len));

   if ~isempty(opt.segType)
      % Check all connections leaving tube

      outFromTube = find(opt.segType == 3);

      for i = 1 : numel(outFromTube)
         segNo = outFromTube(i);

         % Find the two neighbor nodes
         neighborNodes = topo(segNo, :);
         d = 0;
         n = 0;

         for j = 1:2
            % Connected to two nodes - check which segments connected to
            % these two nodes are the tubing, and compute average length
            % over these segments to assign the valve a fraction of the
            % tubing.
            neigh = neighborNodes(j);
            localSegments = find(any(s.topo == neigh, 2));

            for k = 1 : numel(localSegments)
               lseg = localSegments(k);

               if opt.segType(lseg) == 1
                  d = d + len(lseg);
                  n = n + 1;
               end
            end
         end

         assert (n > 0)

         dist(segNo) = d ./ n;
      end
   end
end
