function [trap_ixs, vols] = downstreamTraps(Gt, poro, tstruct, trap_ix)
%
% Return index of all traps downstream of 'trap_ix'.  Also include
% 'trap_ix'.  If output argument 'vols' is provided, the trap volumes will
% also be computed.
%
% SYNOPSIS:
%   function [trap_ixs, vols] = downstreamTraps(Gt, poro, tstruct, trap_ix)
%
% PARAMETERS:
%   Gt      - Top surface grid in which the trap lies
%   poro    - Porosity field (one value per cell in Gt)
%   tstruct - trap structure og Gt, as computed by 'trapAnalysis'
%   trap_ix - index of starting trap
%
% RETURNS:
%   trap_ixs - Indices of all traps downstream of 'trap_ix', as well as
%              'trap_ix' itself
%   vols     - (optional) pore volumes of each downstream trap (one value per
%              entry in 'trap_ixs')
%
% SEE ALSO:
% trapAnalysis

  trap_ixs = [];
  new_traps = trap_ix;
  adj = sparse(eye(size(tstruct.trap_adj,1))); % Identity matrix of same size
                                               % as provided adjacency matrix
  
  while ~isempty(new_traps)
      trap_ixs = [trap_ixs, new_traps]; %#ok MATLAB has no suitable datastructure
      adj = adj * tstruct.trap_adj;
      new_traps = find(adj(trap_ix, :));
  end
 
  % we want to preserve the order in which the traps were encountered, so we
  % call unique with the 'stable' option.
  trap_ixs = unique(trap_ixs, 'stable');
  
  if nargout > 1
      num_traps = numel(trap_ixs);
      vols = zeros(1, num_traps);
      
      for i = 1:num_traps
          vols(i) = computeTrapVolume(Gt, tstruct, poro, trap_ixs(i));
      end
  end
end

  

  
