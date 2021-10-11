function bc = addBCVEM(bc, f, t, g)
%Add boundary condition to (new or existing) BC using function handels,
%to be used in patch-testing and convergence tests in incompVEM.
%
% Syntax is similar to that of addBC, but with a function handle in stead
% of scalar value for each boundary. Each call to addBCVEM adds a set of
% faces with a given type and function handle bc.
%
% SYNOPSIS:
%   bc = addBC(bc, faces, type, function)
%
% PARAMETERS:
%   bc       - Boundary condition structure from a prior call to 'addBCVEM'
%              which will be updated on output or an empty array (bc==[])
%              in which case a new boundary condition structure is created.
%
%   faces    - Global faces in external model for which this boundary
%              condition should be applied.
%
%   type     - Type of boundary condition. Supported values are 'pressure'
%              and 'flux', or cell array of such strings.
%
%   function - Boundary condition function. Interpreted as a pressure p
%              (in units of 'Pa') when type=='pressure' and as a
%              -K \nabla p \cdot n when type=='flux'. One scalar value for
%              each face in 'faces'.
%
%              Note: type=='flux' is only supported in 2D, and the function
%              is interpreted as K \nabla p \cdot n.
%
% SEE ALSO:
%   addBC

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

if isempty(f),
   warning('MRST:addBC', 'Empty list of boundary faces.');
   return;
end

if isempty(bc),
   bc = struct('face', {{}}, ...
               'type', {{}}, ...
               'func', {{}}, ...
               'sat' , []       );
end

%   Validate boundary condition type.
assert(ischar(t), 'Boundary condition type should be a string');
t = lower(t);

assert (all(strcmpi(t, 'pressure') | strcmpi(t, 'flux')), ...
      'Boundary condition type should be either ''pressure'' of ''flux''');

%   Verify that boundary condition is not already set
bc_given = cellfun(@(face) nnz(ismember(face, f)), bc.face);
assert(nnz(bc_given) == 0, ...
  'New boundary condition overlaps with conditions already in bc-struct.');
 
%   Update  bc struct.
nn = numel(bc.face)+1;
bc.face{nn} = f;
bc.type{nn} = t;
bc.func{nn} = g;

end