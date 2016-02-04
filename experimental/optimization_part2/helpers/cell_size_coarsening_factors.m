function [Gt, rock2D, coarsening] = cell_size_coarsening_factors( fmName, cell_size, coarsening_limit )
% Determine the best corasening factors to use for CO2 Atlas datasets, such
% that either a desired number of cells or a desired cell resolution is
% achieved.

    % Huginfmeast is the smallest formation in terms of surface area, with
    % 2264 cells at dx=dy=1000 meters (un-coarsened).
    
    % Find the maximum coarsening factor possible
    c = coarsening_limit;
    
    max_coarsening = [];
    while isempty(max_coarsening)
        try
            [~, rock2D] = getFormationTopGrid( fmName, c );
            % if successful, leave while loop
            fprintf('Max coarsening level possible is %4.0f.\n', c)
            max_coarsening = c;
        catch
            % reduce c if it was too high
            c = c - 1;
        end
    end
    
    % Now refine until cell size constraint is met
    clear rock2D
    c = max_coarsening;
    while c > 0
        [Gt, rock2D] = getFormationTopGrid( fmName, c );
        if sqrt(mean(Gt.cells.volumes(:))) > cell_size+1 % @@
            c = c - 1;
            
            % extra check for some bad combinations
            if strcmpi(fmName,'Sognefjordfm') && c == 19
                c = 18;
            elseif strcmpi(fmName,'Fensfjordfm') && c == 19
                c = 18;
            elseif strcmpi(fmName,'Huginfmeast') && c == 10
                c = 9;
            end
            
        else
            break
        end
    end
    if c == 0
       warning('Un-coarsened grid contains cell sizes larger than your requested size.') 
       c = 1;
    end
    
    fprintf('Coarsening of %2.0f gives %6.0f cells a %6.0f meter cell size.\n', ...
        c, Gt.cells.num, sqrt(mean(Gt.cells.volumes(:))))
    coarsening = c;
end

