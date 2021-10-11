function varargout = mixedSymm(B, C, D, f, g, h, Do, varargin)
%Solve symmetric system of linear eqns using reduction to a mixed system.
%
% SYNOPSIS:
%   [faceFlux, press, lam_n]    = mixedSymm(B, C, D, f, g, h, Do)
%   [faceFlux, press, lam_n]    = mixedSymm(B, C, D, f, g, h, Do, ...
%                                           'pn', pv, ...)
%   [faceFlux, press, lam_n, A] = mixedSymm(...)
%
% DESCRIPTION:
%   Solves the symmetric, hybrid (block) system of linear equations
%
%       [  B   C   D  ] [  v ]     [ f ]
%       [  C'  0   0  ] [ -p ]  =  [ g ]
%       [  D'  0   0  ] [ cp ]     [ h ]
%
%   with respect to the face flux 'vf', the cell pressure 'p', and the face
%   (contact) pressure 'cp' for Neumann boundary.  The system is solved
%   using a reduction to a symmetric mixed system of linear equations
%
%       [  Bm   Cm   Dm ] [ vf   ]     [  fm ]
%       [  Cm'  0    0  ] [ -p   ]  =  [  g  ]
%       [  Dm'  0    0  ] [ cp_n ]     [ h_n ]
%
%   where Bm, Cm are 'the' mixed block matrices, and Dm enforces Neumann
%   BCs (size(Dm) = m x mb, where n is the number of faces and mb is the
%   number of boundary faces).
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
%   A        - Mixed system matrix.
%              OPTIONAL.  Only returned if specifically requested.
%
% SEE ALSO:
%   `tpfSymm`, `schurComplementSymm`, `solveIncompFlow`, `solveIncompFlowMS`.

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

   assert (numel(f) == size(B,1));

   numC   = size(C,  2);
   numF   = size(Do, 2);
   numNF  = size(D,  2);

   %-----------------------------------------------------------------------
   % Compute component matrices. ------------------------------------------
   %
   if ~opt.MixedB,
      Bm = Do.' * B * Do;
      fm = Do.' * f;
   else
      Bm = B;
      fm = f;
   end
   Cm = Do.' * C;
   Dm = Do.' * D;

   %-----------------------------------------------------------------------
   % Form (reduced) mixed system. -----------------------------------------
   %
   %   A x = rhs
   %
   A   = [Bm ,               Cm   , Dm    ; ...
          Cm', sparse(numC,  numC + numNF); ...
          Dm', sparse(numNF, numC + numNF)];

   rhs = [fm; g; h];

   %-----------------------------------------------------------------------
   % Solve system. --------------------------------------------------------
   %
   if opt.Regularize,
      A(numF + 1, numF + 1) = 1;
   end
   x = opt.LinSolve(A, rhs);

   %-----------------------------------------------------------------------
   % Extract output from solution vector. ---------------------------------
   %
   varargout{1} =  x(              1 : numF       );  % faceFlux
   varargout{2} = -x(       numF + 1 : numF + numC);  % press
   varargout{3} =  x(end - numNF + 1 : end        );  % lam_n

   if nargout > 3, varargout{4} = A; end
end
