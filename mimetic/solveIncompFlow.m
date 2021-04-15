function varargout = solveIncompFlow(varargin)
%Solve incompressible flow problem (fluxes/pressures). [DEPRECATED]
%
% SYNOPSIS:
%   state = solveIncompFlow(state, G, S, fluid)
%   state = solveIncompFlow(state, G, S, fluid, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function assembles and solves a (block) system of linear equations
%   defining interface fluxes and cell and interface pressures at the next
%   time step in a sequential splitting scheme for the reservoir simulation
%   problem defined by Darcy's law and a given set of external influences
%   (wells, sources, and boundary conditions).
%
% REQUIRED PARAMETERS:
%   state  - Reservoir and well solution structure either properly
%            initialized from function 'initState', or the results from a
%            previous call to function 'solveIncompFlow' and, possibly, a
%            transport solver such as function 'explicitTransport'.
%
%   G, S   - Grid and (mimetic) linear system data structures as defined by
%            function 'computeMimeticIP'.
%
%   fluid  - Fluid data structure as described by 'fluid_structure'.
%
% OPTIONAL PARAMETERS:
%   wells  - Well structure as defined by function 'addWell'.  May be empty
%            (i.e., W = []) which is interpreted as a model without any
%            wells.
%
%   bc     - Boundary condition structure as defined by function 'addBC'.
%            This structure accounts for all external boundary conditions
%            to the reservoir flow.  May be empty (i.e., bc = []) which is
%            interpreted as all external no-flow (homogeneous Neumann)
%            conditions.
%
%   src    - Explicit source contributions as defined by function
%            'addSource'.  May be empty (i.e., src = []) which is
%            interpreted as a reservoir model without explicit sources.
%
%   rhs    - Supply system right-hand side 'b' directly.  Overrides
%            internally constructed system right-hand side.  Must be a
%            three-element cell array, the elements of which are correctly
%            sized for the block system component to be replaced.
%
%            NOTE: This is a special purpose option for use by code which
%            needs to modify the system of linear equations directly, e.g.,
%            the 'adjoint' code.
%
%   Solver - Which solver mode function 'solveIncompFlow' should employ in
%            assembling and solving the block system of linear equations.
%            String.  Default value: Solver = 'hybrid'.
%
%            Supported values are:
%              - 'hybrid' --
%                   Assemble and solve hybrid system for interface
%                   pressures.  System is eventually solved by Schur
%                   complement reduction and back substitution.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','hybrid') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
%              - 'mixed' --
%                   Assemble and solve a hybrid system for interface
%                   pressures, cell pressures and interface fluxes. System
%                   is eventually reduced to a mixed system as per function
%                   'mixedSymm'.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','mixed') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
%              - 'tpfa' --
%                   Assemble and solve a cell-centred system for cell
%                   pressures.  Interface fluxes recovered through back
%                   substitution.
%
%                   The system 'S' must in this case be assembled by
%                   passing option pair ('Type','mixed') or option pair
%                   ('Type','comp_hybrid') to function 'computeMimeticIP'.
%
%   LinSolve -
%            Handle to linear system solver software to which the fully
%            assembled system of linear equations will be passed.  Assumed
%            to support the syntax
%
%                        x = LinSolve(A, b)
%
%            in order to solve a system Ax=b of linear equations.
%            Default value: LinSolve = @mldivide (backslash).
%
%   MatrixOutput -
%            Whether or not to return the final system matrix 'A' to the
%            caller of function 'solveIncompFlow'.
%            Logical.  Default value: MatrixOutput = FALSE.
%
% RETURNS:
%   state - Update reservoir and well solution structure with new values
%           for the fields:
%              - pressure -- Pressure values for all cells in the
%                            discretised reservoir model, 'G'.
%              - facePressure --
%                            Pressure values for all interfaces in the
%                            discretised reservoir model, 'G'.
%              - flux     -- Flux across global interfaces corresponding to
%                            the rows of 'G.faces.neighbors'.
%              - A        -- System matrix.  Only returned if specifically
%                            requested by setting option 'MatrixOutput'.
%
%              - wellSol  -- Well solution structure array, one element for
%                            each well in the model, with new values for
%                            the fields:
%                              - flux     -- Perforation fluxes through all
%                                            perforations for corresponding
%                                            well.  The fluxes are
%                                            interpreted as injection
%                                            fluxes, meaning positive
%                                            values correspond to injection
%                                            into reservoir while negative
%                                            values mean
%                                            production/extraction out of
%                                            reservoir.
%                              - pressure -- Well bottom-hole pressure.
%
% NOTE:
%
%   This function is deprecated and will be removed in a future release of MRST
%   Please switch to replacement function 'incompMimetic'.
%
%   If there are no external influences, i.e., if all of the structures
%   'W', 'bc', and 'src' are empty and there are no effects of gravity, and
%   no system right hand side has been supplied externally, then the input
%   state is returned unchanged and a warning is printed in the command
%   window.  This warning is printed with message ID
%
%           'incompMimetic:DrivingForce:Missing'
%
% SEE ALSO:
%   `computeMimeticIP`, `addBC`, `addSource`, `addWell`, `initSimpleFluid`
%   `initState`, `solveIncompFlowMS`, `incompMimetic`.

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

   persistent CALLED;

   if isempty(CALLED), CALLED = false; end

   if ~CALLED,
      warning('solveIncompFlow:Deprecated', ...
             ['Function ''%s'' is deprecated and will be ', ...
              'removed in a future release of MRST\n', ...
              'Please switch to function ''incompMimetic''.'], ...
              mfilename);

      CALLED = true;
   end

   [varargout{1:nargout}] = incompMimetic(varargin{:});
end
