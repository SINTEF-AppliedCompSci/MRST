function table = printTrapInfoTable(Gt, z_spill, traps, ...
                                    sort_column, max_lines, rock2D)
% Print basic information about the traps.
% 
% SYNOPSIS:
%   printTrapInfoTable(Gt, z_spill, traps)
%
% PARAMETERS:
%   Gt          - 2D top-level grid with computed geometry
%   z_spill     - z spill height for each trap (vector)
%   traps       - for each cell, index of trap it belongs to (or zero if not
%                 in any trap)
%   sort_column - Which column to sort the table by:
%                 1: index of trap
%                 2: number of cells in trap
%                 3: volume of trap
%                 Default is by ascending order.  A negative sign indicates
%                 to sort in descending order.
%   max_lines   - If nonzero, specifies the number of lines of the table to
%                 print (rest will be truncated)
%   rock2D      - Rock properties corresponding to the grid Gt.  Only
%                 porosity values are used

% Traps are numbered 1 - number of traps, so the highest index indicates how
% many traps there are in total.

    num_traps = max(traps)
    fprintf('There are %d topologically distinct traps.\n', max(traps));

    % Computing information to go into table
    table = compute_table(Gt, z_spill, traps, rock2D);
    table = sortrows(table, sort_column);
    total_volume = sum(table(:,3));

    % Printing table
    if (max_lines == 0)
        lines_to_print = num_traps;
    else
        lines_to_print = min(num_traps, max_lines);
    end
    fprintf('\nTrap Ix Cells Levels  Volume (%% of total)\n\n')
    
    for i=1:lines_to_print
        fprintf('%7d%6d%10.2e (%3.1f)\n', ...
                [table(i, :) table(i,3)/total_volume * 100]);
    end

    % Printing total volume
    fprintf('\n\nTotal trap volume: d%10.2e\n\n', sum(table(:,3)));

end % _______ end printTrapInfoTable ______
    
%_______________________________________________________________________
function res = compute_table(Gt, spill_z, traps, rock2D)
    
    num_traps = max(traps);
    res = zeros(num_traps, 3);
    for i = 1:num_traps
        covered_cells = find(traps==i);
        res(i,:) = [i ...                      % trap index
                    size(covered_cells,1)  ... % number of cells in trap
                    compute_volume(Gt, spill_z(i), covered_cells, rock2D)];
    end
end

%_______________________________________________________________________
function vol = compute_volume(Gt, z_spill_value, covered_cells, rock2D)
    cell_heights =  -Gt.cells.z(covered_cells) + z_spill_value;
    vol = Gt.cells.volumes(covered_cells)' * (cell_heights .* ...
                                              rock2D.poro(covered_cells));
    
end

%% DEPRECATED CODE BELOW
%_______________________________________________________________________
% function res = max_cell_level(trap_level)
% % Find the highest traplevel for each cell
%     num_cells = size(trap_level{1}, 1);
%     max_level = size(trap_level, 2) - 1;
%     res = zeros(num_cells, 1); % We must assume there's at least one cell in
%                                % the cell array
%     for i = 1:max_level
%         % cells who are assigned a trap at this level should be seen as
%         % belonging to this level (or higher)
%         res(find(sum(trap_level{i},2)>0)) = i;
%     end
% end
