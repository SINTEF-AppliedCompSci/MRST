function n = translateOptionsAMGCL(name, value)
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
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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

    if ~iscellstr({name, value})
       error('OptionArgType:Unsupported', ...
            ['Inputs ''name'' and ''value'' must ', ...
             'both be character vectors']);
    end

    switch lower(name)
        case 'preconditioner'
            n = translate_precondtioner(value);

        case 'coarsening'
            n = translate_coarsening_scheme(value);

        case 'relaxation'
            n = translate_relaxation_scheme(value);

        case 'solver'
            n = translate_solver_name(value);

        otherwise
            error('Option:Unknown', 'Unknown option: ''%s''', name);
    end
end

%--------------------------------------------------------------------------

function n = translate_precondtioner(value)
   switch lower(value)
      case 'amg',        n = 1;
      case 'relaxation', n = 2;
      case 'dummy',      n = 3;

      otherwise
         error('Unknown preconditioner option: ''%s''', value);
   end
end

%--------------------------------------------------------------------------

function n = translate_coarsening_scheme(value)
   switch lower(value)
      case 'smoothed_aggregation', n = 1;
      case 'ruge_stuben',          n = 2;
      case 'aggregation',          n = 3;
      case 'smoothed_aggr_emin',   n = 4;

      otherwise
         error('Unknown coarsening option: ''%s''', value)
   end
end

%--------------------------------------------------------------------------

function n = translate_relaxation_scheme(value)
   switch lower(value)
      case 'spai0',         n = 1;
      case 'gauss_seidel',  n = 2;
      case 'ilu0',          n = 3;
      case 'iluk',          n = 4;
      case 'ilut',          n = 5;
      case 'damped_jacobi', n = 6;
      case 'spai1',         n = 7;
      case 'chebyshev',     n = 8;

      otherwise
         error('Unknown relaxation option: ''%s''', value);
   end
end

%--------------------------------------------------------------------------

function n = translate_solver_name(value)
   switch lower(value)
      case 'bicgstab',  n = 1;
      case 'cg',        n = 2;
      case 'bicgstabl', n = 3;
      case 'gmres',     n = 4;
      case 'lgmres',    n = 5;
      case 'fgmres',    n = 6;
      case 'idrs',      n = 7;

      otherwise
         error('Unknown:Solver', 'Unknown solver name: ''%s''', value);
   end
end
