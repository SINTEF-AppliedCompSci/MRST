function [n, name, description, parameter] = translateOptionsAMGCL(name, value)
%Translate AMGCL String Option Value to Integer Code
%
% SYNOPSIS:
%    n = translateOptionsAMGCL(name, value)
%
% PARAMETERS:
%   name  - Option name.  Character vector.  Must be one of 'precondtioner',
%           'coarsening', 'relaxation', or 'solver'.
%
%   value - Option value.  String (character vector) pertaining to 'name'.
%
% RETURNS:
%   n - Integer option value known to underlying MEX gateway.
%
% SEE ALSO:
%   `callAMGCL`, `amgcl_matlab`, `getAMGCLMexStruct`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
    assert(ischar(name), 'First input must be a name');

    switch lower(name)
        case 'preconditioner'
            choices = {'amg', 'relaxation', 'dummy'};
            descriptions = {'Algebraic multigrid', ...
                            'Relaxation as preconditioner', ...
                            'Dummy (do nothing)'};
            parameters = cell(1, 3);

        case 'coarsening'
            choices = {'smoothed_aggregation', 'ruge_stuben', ...
                       'aggregation', 'smoothed_aggr_emin'};
            descriptions = {'Smoothed aggregation', ...
                            'Ruge-Stuben / classic AMG coarsening', ...
                            'Aggregation with constant interpolation', ...
                            'Smoothed aggregation (energy minimizing)'};
            psg = {'aggr_eps_strong', 'aggr_over_interp', 'aggr_relax'};
            prs = {'rs_eps_strong', 'rs_trunc', 'rs_eps_trunc'};
            pag = {'aggr_eps_strong', 'aggr_over_interp', 'aggr_relax'};
            pem = pag;
            parameters = {psg, prs, pag, pem};

        case 'relaxation'
            choices = {'spai0', 'gauss_seidel', 'ilu0', 'iluk',...
                       'ilut', 'damped_jacobi', 'spai1', 'chebyshev'};
            descriptions = {'Sparse approximate inverse of order 0', ...
                            'Gauss-Seidel smoothing.', ...
                            'Incomplete LU-factorization with zero fill-in - ILU(0)', ...
                            'Incomplete LU-factorization of order k - ILU(k)', ...
                            'Incomplete LU-factorization with thresholding - ILU(t)', ...
                            'Damped Jacobi smoothing', ...
                            'Sparse approximate inverse of order 1', ...
                            'Chebyshev smoothing'};
            plu = {'ilu_damping'};
            plk = {'ilu_damping, ilu_k parameter'};
            plt = {'ilu_damping, ilut_tau'};
            pj =  {'jacobi_damping'};
            pch = {'chebyshev_degree', 'chebyshev_lower', 'chebyshev_power_its'};
            parameters = {{}, {}, plu, plk, plt, pj, {}, pch};
        case 'solver'
            choices = {'bicgstab', 'cg', 'bicgstabl', 'gmres', 'lgmres', 'fgmres', 'idrs'};
            descriptions = {'Biconjugate gradient stabilized method.', ...
                            'Conjugate gradient method.', ...
                            'Biconjugate gradient stabilized method (l variant)', ...
                            'Generalized minimal residual method', ...
                            'Generalized minimal residual method', ...
                            'Generalized minimal residual method (f variant)', ...
                            'Induced Dimension Reduction method' ...
                            };
            pbcg = {'bicgstabl_convex', 'bicgstabl_delta', 'bicgstabl_l'};
            pgm =  {'gmres_m'};
            plgm = {'lgmres_always_reset', 'lgmres_k', 'lgmres_store_av'};
            pdrs = {'idrs_omega', 'idrs_replacement', 'idrs_s'};
            parameters = {{}, {}, pbcg, pgm, plgm, {}, pdrs};
        otherwise
            error('Option:Unknown', 'Unknown option: ''%s''', name);
    end
    [n, name, description, parameter] = translate_option(name, choices, descriptions, parameters, value);
end


function [index, name, description, option] = translate_option(category, names, descriptions, options, value)
   index = mapArgument(value, names, category);
   if ischar(index)
       name = names;
       description = descriptions;
       option = options;
   else
       name = names{index};
       description = descriptions{index};
       option = options{index};
   end
end

function index = mapArgument(value, types, groupName)
   if isempty(value)
       % Get everything
       index = ':';
   elseif ischar(value)
       % String corresponding to valid option
       index = find(strcmpi(types, value));
       if isempty(index)
           error(['Unknown:', groupName], 'Unknown %s name: ''%s''', ...
                                            groupName, value);
       end
   else
       % Numeric index
       index = value;
       assert(isnumeric(value));
       assert(isnumeric(value) <= numel(types));
   end
end
