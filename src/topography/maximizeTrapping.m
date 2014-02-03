function [trees, v] = maximizeTrapping(G, varargin)
%Find the N best injection trees, with optional compensation for overlap
%
% SYNOPSIS:
%   
%
% DESCRIPTION:
%   trees = maximizeTrapping(G)
%   [trees, volumes] = maximizeTrapping(G, 'n', 5)
%
% REQUIRED PARAMETERS:
%   G   - Top surface grid
%
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   N   - The number of trapping trees to consider. This will, for a
%         most purposes, be the number of injection wells to be drilled for
%         a CO2 migration study.
%
%   res - Output from trapAnalysis. Will be generated if not supplied.
% 
%   calculateAll - Should all traps be considered as starting points? If
%                  not, leaf nodes will be the starting point of the
%                  algorithm.
%
%   removeOverlap - If enabled, traps belonging to tree number 1, ..., m-1
%                   will be set to zero volume when calculating tree number
%                   m.
%
% RETURNS:
%   trees     - N by 1 struct array of trees structs, with fields
%                   - root: Trap index of the tree root.
%                   - traps: Array of traps downstream from root (including
%                            the root trap itself)
%                   - value: Total value of the tree. If N >= number of
%                   traps and removeOverlap is enabled, this will sum to
%                   total trap volume over all trees.

%{
#COPYRIGHT#
%}

    opt = struct('res',             [],...
                 'n',               inf,...
                 'removeOverlap',   true, ...
                 'calculateAll',    false);
 
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.res)
        opt.res = trapAnalysis(G, true);
    end
    
    v = trapVolumes(G, opt.res, 1:max(opt.res.traps));

    % All root nodes - traps that are downstream from every other trap
    if opt.calculateAll
        rn = 1:max(opt.res.traps);
    else
        rn = rootNotes(opt.res.trap_adj');
    end
    Nt = min(opt.n,numel(rn));
    
    trees = struct('root',  [], ...
                   'traps', [], ...
                   'value', cellfun(@(x) 0, cell(1,Nt), 'UniformOutput', false));
    % Initialize
    i = 1;
    
    downstream = arrayfun( @(ind) find(dfs(opt.res.trap_adj, ind) > -1), rn, ...
                                                  'UniformOutput', false);
    rn_vols = cellfun(@(x) sum(v(x)), downstream);
    done = false(size(rn_vols));
    while i <= Nt
        [vol, index] = max(rn_vols);
        
        ds = downstream{index};
        trees(i).value = vol;
        trees(i).traps = ds;
        trees(i).root = rn(index);
        if opt.removeOverlap
            v(ds) = 0;
        end
        % Flag node as done - we will not use it as a new root node 
        done(index) = true;
        
        rn_vols = cellfun(@(x) sum(v(x)), downstream);
        rn_vols(done) = 0;
        i = i + 1;
    end
end



function c = rootNotes(A)
    c = find(sum(A,2) == 0);
end
