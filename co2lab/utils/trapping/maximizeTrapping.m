function [trees, vols] = maximizeTrapping(G, varargin)
%Find the N best injection trees, with optional compensation for overlap
%
% SYNOPSIS:
%   trees            = maximizeTrapping(G)
%   [trees, volumes] = maximizeTrapping(G, 'n', 5)
%
% REQUIRED PARAMETERS:
%   G   - Top surface grid
%
%
% OPTIONAL PARAMETERS:
%   N   - The number of trapping trees to consider. This will, for a
%         most purposes, be the number of injection wells to be drilled for
%         a CO2 migration study.
%
%   res - Output from trapAnalysis. Will be generated if not supplied.
% 
%   calculateAll - Should all traps be considered as starting points? If
%         not, leaf nodes will be the starting point of the algorithm.
%
%   removeOverlap - If enabled, traps belonging to tree number 1, ..., m-1
%         will be set to zero volume when calculating tree number m.
%
% RETURNS:
%   trees - N by 1 struct array of trees structs, with fields
%           - root:  Trap index of the tree root.
%           - traps: Array of traps downstream from root (including
%                    the root trap itself)
%           - value: Total value of the tree. If N >= number of
%                    traps and removeOverlap is enabled, this will sum to
%                    total trap volume over all trees.
%   vols  - Individual trap volumes
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

    opt = struct('res',             [],...
                 'n',               inf,...
                 'removeOverlap',   true, ...
                 'calculateAll',    false);
 
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.res)
        opt.res = trapAnalysis(G, true);
    end
    
    v = volumesOfTraps(G, opt.res, 1:max(opt.res.traps));
    vols = v;
    
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
    
    % downstream = arrayfun( @(ind) find(dfs(opt.res.trap_adj, ind) > -1), rn, ...
    %                                               'UniformOutput', false);
    downstream = establishDownstreamTraps(opt.res.trap_adj, rn);
    
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
% ============================================================================

function dstr = establishDownstreamTraps(A, start_ixs)
% Returns cell array with one entry per trap indexed in 'start_ixs', listing
% the chain of traps it spills into (including itself)
    num_traps = size(A,1); % also equal to size(A,2)
    M  = spones(A + speye(num_traps));
    M1 = M;
    M2 = spones(M1 * M);
    while ~isequal(M1, M2)
        M1 = M2;
        M2 = spones(M1 * M);
    end
    dstr = arrayfun(@(x) find(M2(x,:))', start_ixs, 'UniformOutput', false);
end

% ----------------------------------------------------------------------------
function c = rootNotes(A)
    c = find(sum(A,2) == 0);
end
