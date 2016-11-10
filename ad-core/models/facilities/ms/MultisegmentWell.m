classdef MultisegmentWell < SimpleWell
    properties
        operators
    end
    
    methods
        function well = MultisegmentWell(W, varargin)
            well = well@SimpleWell(W);
            assert(W.isMS)
            if nargin > 1
                well = merge_options(well, varargin{:});
            end
            [nn, ns] = deal(numel(w.nodes.depth), numel(w.segments.length));
            C = sparse((1:ns)'*[1 1], w.segments.topo(1:end,:), ones(ns,1)*[1, -1], ns, nn);
            well.operators.grad = @(x)-C*x;
            well.operators.div  = @(x)C'*x;
            well.operators.segmentUpstr = @(flag, val)segmentUpstreamValue(flag, val, w.segments.topo);
            well.operators.C = C;
            aver = sparse((1:ns)'*[1 1], w.segments.topo(1:end,:), ones(ns,1)*[1, ...
                                1], ns, nn);
            well.operators.aver = bsxfun(@rdivide, aver, sum(aver, 2));
        end
    end
end

function v = segmentUpstreamValue(flag, val, topo)
    ix = flag.*topo(:,1) + ~flag.*topo(:,2);
    if ix(1) == 0
        ix(1) = topo(1,2);
    end
    v = val(max(ix-1,1));
end
