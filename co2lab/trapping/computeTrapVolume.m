function vol = computeTrapVolume(Gt, tstruct, poro, trap_ixs)
% Compute total pore volume of traps
%
% SYNOPSIS:
%   vol = computeTrapVolume(Gt, tstruct, poro, trap_ix)
%
% PARAMETERS:
%   Gt       - top surface grid
%   tstruct  - trap structure of 'Gt', as produced by 'trapAnalysis'
%   poro     - porosity (one value per cell in Gt)
%   trap_ixs - indices of traps for which volume should be computed
%
% RETURNS:
%   vol - total pore volumes of the requested traps
%
% SEE ALSO:
% trapAnalysis

  vol = zeros(numel(trap_ixs),1);
    
  for i = 1:numel(trap_ixs)
      % find index of all trap cells in this trap
      cix = find(tstruct.traps == trap_ixs(i));
      
      trap_heights = tstruct.trap_z(trap_ixs(i)) - Gt.cells.z(cix);
      trap_heights = min(trap_heights, Gt.cells.H(cix));
      
      cell_vols = trap_heights .* Gt.cells.volumes(cix) .* poro(cix);
      vol(i) = sum(cell_vols);
  end
end

  