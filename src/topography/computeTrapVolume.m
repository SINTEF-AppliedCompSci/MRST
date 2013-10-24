function vol = computeTrapVolume(Gt, tstruct, poro, trap_ixs)
% Compute total pore volume of the traps whose indices are given in 'trap_ixs'
% in the top surface grid Gt.
% 
% SYNOPSIS:
%   function vol = computeTrapVolume(Gt, tstruct, poro, trap_ix)
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
      cell_vols = (tstruct.trap_z(trap_ixs(i)) - Gt.cells.z(cix)) ...
          .* Gt.cells.volumes(cix) ...
          .* poro(cix);
      vol(i) = sum(cell_vols);
  end
end

  