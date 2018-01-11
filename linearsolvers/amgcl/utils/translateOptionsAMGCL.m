function n = translateOptionsAMGCL(name, value)
    switch lower(name)
        case 'preconditioner'
            switch(value)
                case 'amg'
                    n = 1;
                case 'relaxation'
                    n = 2;
                case 'dummy'
                    n = 3;
                otherwise
                    error('Unknown preconditioner option.')
            end
        case 'coarsening'
            switch(value)
                case 'smoothed_aggregation'
                    n = 1;
                case 'ruge_stuben'
                    n = 2;
                case 'aggregation'
                    n = 3;
                case 'smoothed_aggr_emin'
                    n = 4;
                otherwise
                    error('Unknown coarsening option.')
            end
        case 'relaxation'
            switch(value)
                case 'spai0'
                    n = 1;
                case 'gauss_seidel'
                    n = 2;
                case 'ilu0'
                    n = 3;
                case 'iluk'
                    n = 4;
                case 'ilut'
                    n = 5;
                case 'damped_jacobi'
                    n = 6;
                case 'spai1'
                    n = 7;
                case 'chebyshev'
                    n = 8;
                otherwise
                    error('Unknown relaxation option.')
            end
        case 'solver'
            switch(value)
                case 'bicgstab'
                    n = 1;
                case 'cg'
                    n = 2;
                case 'bicgstabl'
                    n = 3;
                case 'gmres'
                    n = 4;
                case 'lgmres'
                    n = 5;
                case 'fgmres'
                    n = 6;
                case 'idrs'
                    n = 7;
                otherwise
                    error('Unknown solver option.')
            end
        otherwise
            error('Unknown option');
    end
end