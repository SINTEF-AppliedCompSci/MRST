function ustream_trap_indices = identify_upstream_traps(start_traps, connectivity)
%
% SYNOPSIS:
%   function ustream_trap_indices = identify_upstream_traps(start_traps, connectivity)
%
% PARAMETERS:
%   start_traps  - index(ices) of trap(s) for which we will identify the upstream traps
%   connectivity - adjacency matrix for traps
%
% RETURNS:
%   ustream_trap_indices - indices of start traps and all upstream traps, in the order
%                          encountered (duplicates removed). 

encountered_traps = start_traps; % start
current_traps = start_traps;

% locating all upstream traps
while current_traps
    current_traps = find(sum(connectivity(current_traps, :),1));
    encountered_traps= [encountered_traps, current_traps];
end

% removing duplicates, if any (but preserve relative order)

[~, ix] = unique(encountered_traps, 'first');
ustream_trap_indices = encountered_traps(sort(ix));

