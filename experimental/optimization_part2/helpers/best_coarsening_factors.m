function [Gt, rock2D, coarsening] = best_coarsening_factors( fmName, varargin )
% Determine the best corasening factors to use for CO2 Atlas datasets, such
% that either a desired number of cells or a desired cell resolution is
% achieved.

    opt.max_cell_size = 3000 * meter;
    opt.cell_size = 3000 * meter;
    %opt.max_num_cells = 10000;
    opt = merge_options(opt, varargin{:});

    % We start with a very coarse resolution, and asses number of cells and
    % cell resolution. We refine until desired contraints are satified.
    c = 9; % the max coarsening which Huginfmeast can handle 
    while c > 0
        [Gt, rock2D] = getFormationTopGrid( fmName, c );
        if sqrt(mean(Gt.cells.volumes(:))) > opt.max_cell_size %|| Gt.cells.num > opt.max_num_cells
            c = c - 1;
        else
            break
        end
    end
    if c == 0
       warning('Un-refined grid does not meet constraints.') 
       c = 1;
    end
    
    fprintf('Coarsening of %2.0f gives %6.0f cells a %6.0f meter cell size.\n', ...
        c, Gt.cells.num, sqrt(mean(Gt.cells.volumes(:))))
    coarsening = c;
end

