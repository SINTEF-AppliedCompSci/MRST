function [fill_heights, sat_field, escaped_vol] = ...
    computeTrappedDistribution(Gt, tstruct, cell_ix, vol, poro)
% 
% Compute the final distribution of an injected quantity of CO2, only taking
% topography/structural trapping into account, and assuming infinitely slow flow.
%
% SYNOPSIS:
%   function [fill_heights, hfield, escaped_vol] =  
%      computeTrappedDistribution(Gt, tstruct, cell_ix, vol, poro)
%
% PARAMETERS:
%   Gt      - top surface grid
%   tstruct - trapping structure for 'Gt', as computed by 'trapAnalysis'
%   cell_ix - index of the cell into which CO2 is injected
%   vol     - volume of injected CO2
%   poro    - porosity of each cell in 'Gt'.
%
% RETURNS:
%   fill_heights - one value per trap in the grid.  Expressing to which
%                  z-height each trap has been filled.  A zero value means
%                  that this trap has gone untouched.
%   sat_field    - saturation field (one value per cell)
%   escaped_vol  - the amount of injected volume that has not been trapped
%                  (i.e. that has 'escaped')
%
% EXAMPLE:
%
% SEE ALSO:
% trapAnalysis
% 
  num_traps = numel(tstruct.trap_z);
  fill_heights = zeros(num_traps, 1); % initially, no trap has received any CO2

  escaped_vol = vol; % Initially assume that all injected volume will escape
  sat_field = zeros(Gt.cells.num, 1);
  
  start_trap = tstruct.trap_regions(cell_ix);
  if start_trap == 0 
      % cell does not belong to any interior trap spill region.  All injected
      % fluid escapes
      return;
  end
  
  % Computing the traps involved and their degree of fill
  [fill_heights, escaped_vol] = ...
      distribute_downwards(Gt, tstruct, start_trap, vol, poro, fill_heights);
  
  % Finally, computing the saturation field
  for t_ix = find(fill_heights > 0)'
      h = fill_heights(t_ix);
      cixs = find(tstruct.traps == t_ix);
      sat = Gt.cells.volumes(cixs) .* (h - Gt.cells.z(cixs)) ./ Gt.cells.H(cixs);
      sat = max(0, sat); % eliminates cells positioned higher than 'h'
      sat_field(cixs) = sat;
  end
  
end

% ==============================================================================
function [fill_heights, escaped] = distribute_downwards(Gt, tstruct, trap_ix, vol, poro, fill_heights)

    % compute how much of the injected volume can be trapped in the present trap
    tvol = empty_trap_volume(Gt, tstruct, trap_ix, poro, fill_heights(trap_ix));
    %fprintf('>> Trap volume of trap %d was %f\n', trap_ix, tvol);

    if tvol >= vol  % this trap can contain all the injected volume
        member_cells = find(tstruct.traps == trap_ix);
        fill_heights(trap_ix) = ...
            compute_exact_fill_height(Gt.cells.volumes(member_cells), ...
                                      Gt.cells.z(member_cells), ...
                                      poro(member_cells), ...
                                      tstruct.trap_z(trap_ix), ...
                                      vol, ...
                                      fill_heights(trap_ix));
        escaped = 0;
        return;
    end
    
    % if we got here, the trap was not large enough to hold all the injected
    % volume.  We must keep track of where the surplus volume ended up.
    fill_heights(trap_ix) = tstruct.trap_z(trap_ix); % ... and completely filled
    
    remaining_vol = vol - tvol;
        
    % determining downstream trap(s)
    dtraps = find(tstruct.trap_adj(trap_ix,:));
    
    if isempty(dtraps) 
        % no downstream trap - remainder escapes
        escaped = remaining_vol;
    else

        escaped = 0; % we will accumulate into this below
        for i = dtraps  

            %  we assume that the remaining volume is divided evenly by all
            %  immediately downstream traps (if more than one...)
            [fill_heights, esc] = ...
                distribute_downwards(Gt, tstruct, i, remaining_vol/numel(dtraps), ...
                                     poro, fill_heights);
            escaped = escaped + esc;
        end
    end
end

% ==============================================================================
function vol = volume_downto_z(areas, top_zvals, bottom_z, poro)

    cell_vols = areas .* (bottom_z - top_zvals) .* poro;
    cell_vols = max(cell_vols, 0); % ignore cells whose z-value is lower than bottom_z
                                   % anyway
    vol = sum(cell_vols);
end


% ==============================================================================
function vol = empty_trap_volume(Gt, tstruct, trap_ix, poro, z_cap)
% 'zcap' > 0 represents a partially filled trap.  The trap is already filled up to
% z-value 'zcap'
    ixs = find(tstruct.traps == trap_ix);
    zvals = max(Gt.cells.z(ixs), z_cap);
    
    vol = volume_downto_z(Gt.cells.volumes(ixs), zvals, ...
                          tstruct.trap_z(trap_ix), poro(ixs));
end

% ==============================================================================
function fheight = compute_exact_fill_height(areas, top_zvals, poro, ubound, ...
                                             volume, zcap)

    % 'zcap' represents the part of the trap already filled (if any)
    %
    % monotonous function that _should_ cross zero if the
    % 'compute_exact_fill_height' function is only invoked when it should be
    % (i.e. when 'volume' is positive but small enough that it doesn't
    % completely fill the remaining volume of the trap.
    
    vol_diff = @(X) volume_downto_z(areas, max(top_zvals, zcap), X, poro) - volume;
    fheight = fzero(vol_diff, 0.5 * (min(top_zvals) + ubound));
    
    
    % Function minimization approach (probably less precise, due to the
    % squaring, but perhaps more robust, due to explicit bounds)
    % obj_fun = @(X) (volume_downto_z(areas, top_zvals, X, poro) - volume)^2;
    % fheight = fminbnd(obj_fun, min(top_zvals), ubound);
    
end
