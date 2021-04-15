function varargout = tpfSymm(B, C, D, f, g, h, Do, varargin)
%Solve symmetric system of linear eqns using reduction to two-point system.
%
% SYNOPSIS:
%   [faceFlux, press, lam_n]    = tpfSymm(B, C, D, f, g, h, Do)
%   [faceFlux, press, lam_n]    = tpfSymm(B, C, D, f, g, h, Do, ...
%                                         'pn', pv, ...)
%   [faceFlux, press, lam_n, A] = tpfSymm(...)
%
% DESCRIPTION:
%   Solves the symmetric, hybrid (block) system of linear equations
%
%       [  B   C   D  ] [  v ]     [ f ]
%       [  C'  0   0  ] [ -p ]  =  [ g ]
%       [  D'  0   0  ] [ cp ]     [ h ]
%
%   with respect to the face flux 'vf', the cell pressure 'p', and the face
%   (contact) pressure 'cp' for neumann boundary.  The system is solved
%   using a reduction to two-point flux system of linear equations
%
%        A [p cp_n] = rhs
%
%   where A is size m-by-m, m = number of cells + number of Neumann faces.
%   If B is not diagonal, off-diagonal terms are ignored.
%
% PARAMETERS:
%   B       - Matrix B  in the block system.  Assumed SPD.
%   C       - Block 'C' of the block system.
%   D       - Block 'D' of the block system, reduced so as to only apply to
%             the Neumann faces (i.e., external faces for which the
%             pressure is not prescribed).
%   f       - Block 'f' of the block system right hand side.
%   g       - Block 'g' of the block system right hand side.
%   h       - Block 'h' of the block system right hand side, reduced so as
%             to apply only to Neumann faces.
%   Do      - Matrix mapping face-fluxes to half-face-fluxes.
%   'pn'/pv - List of name/value pairs describing options to solver.
%             The following options are supported:
%               - 'Regularize' -- TRUE/FALSE, whether or not to enforce
%                                 p(1)==0 (i.e., set pressure zero level).
%                                 This option may be useful when solving
%                                 pure Neumann problems.
%
%               - 'MixedB'     -- TRUE/FALSE (Default FALSE). Must be set
%                                 to TRUE if input 'B' is actually the
%                                 mixed dicsretization matrix Bm (and not
%                                 the hybrid). Useful for multiscale when
%                                 Bm is formed directly.
%               - 'LinSolve'   -- Function handle to a solver for the
%                                 system above.  Assumed to support the
%                                 syntax
%
%                                    x = LinSolve(A, b)
%
%                                 to solve a system Ax=b of linear
%                                 equations.
%                                 Default value: LinSolve = @mldivide
%                                 (backslash).
%
% RETURNS:
%   faceFlux - The face fluxes 'vf'.
%   press    - The cell pressure 'p'.
%   lam_n    - The contact pressure 'cp' on Neumann faces.
%   A        - TPF system matrix.  OPTIONAL.  Only returned if specifically
%              requested.
%
% SEE ALSO:
%   `mixedSymm`, `schurComplementSymm`, `solveIncompFlow`, `solveIncompFlowMS`.

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


   opt = struct('Regularize', false, ...
                'MixedB',     false, ...
                'LinSolve',   @mldivide);
   opt = merge_options(opt, varargin{:});

   numC  = size(C , 2);
   numF  = size(Do, 2);
   numNF = size(D , 2);

   %-----------------------------------------------------------------------
   % Compute component matrices. ------------------------------------------
   %
   if ~opt.MixedB,
      Bm = Do.' * B * Do;
   else
      Bm = B;
   end
   CmDm = Do.' * [C, D];
   fm   = Do.' * f;

   %-----------------------------------------------------------------------
   % Form TPF system. -----------------------------------------------------
   %
   %    A x = rhs
   %
   % Assumes 'Bm' diagonal.
   %
   BiCmDmfm = sparse(1 : numF, 1 : numF, diag(Bm)) \ [CmDm, fm];
   A   = CmDm.' * BiCmDmfm(:,1:end-1);
   rhs = [g; h] - CmDm.'*BiCmDmfm(:,end);

   %-----------------------------------------------------------------------
   % Solve system. --------------------------------------------------------
   %
   if opt.Regularize,
      A(1,1) = 2 * A(1,1);
   end
   x = opt.LinSolve(A, rhs);

   %-----------------------------------------------------------------------
   % Extract output from solution vector. ---------------------------------
   %
   varargout{1} = BiCmDmfm * [x; 1];          % faceFlux
   varargout{2} = x(1 : numC);                % press
   varargout{3} = x(end - numNF + 1 : end);   % lam_n

   if nargout > 3, varargout{4} = A; end
end
