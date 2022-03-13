function W = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, varargin)
%Construct coarse system component matrices for well contributions.
%
% SYNOPSIS:
%   W = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W)
%   W = generateCoarseWellSystem(G, S, CG, CS, mob, rock, W, ...
%                                'pn1', pv1, ...)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   S       - System structure describing the underlying fine grid model,
%             particularly the individual cell flux inner products.
%
%   CG      - Coarse grid structure as defined by function
%             generateCoarseGrid.
%
%   CS      - Coarse system structure as defined by function
%             generateCoarseSystem.
%
%   mob     - Total mobility.  One scalar value for each cell in the
%             underlying (fine) model.
%
%   rock    - Rock data structure.  May be empty, but must contain valid
%             field 'rock.perm' if the synthetic driving source term for
%             basis functions is 'perm'.  Similarly, if the driving source
%             term is 'poros', then the rock data structure must contain a
%             valid field 'rock.poro'.
%
%   W       - Well structure as defined by function addWell.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%                - src -- Explicit source terms in the underlying fine grid
%                         model which must be taken into account when
%                         generating the basis functions.  Note that basis
%                         functions in blocks containing an explicit source
%                         term will be generated based solely on the
%                         explicit source.  The values in 'weight'
%                         pertaining to such blocks will be ignored.
%
%                         Must be a source data structure as defined by
%                         function 'addSource'.  Default value is [] (an
%                         empty array), meaning that no explicit sources
%                         are present in the model.
%
%                - OverlapWell --
%                         Number of fine grid cells by which to extend
%                         support of the well basis function about the well
%                         bore.
%
%                - OverlapBlock --
%                         Number of fine grid cells by which to extend
%                         support of the well basis function about the
%                         reservoir coarse block.
%
% RETURNS:
%   W - Updated well structure having new field W(i).CS defined as
%       W(i).CS.basis : Cell array of SPARSE input vectors from which well
%                       flux basis function matrix may be formed.
%       W(i).CS.basisP: Cell array of SPARSE input vectors from which well
%                       pressure basis function matrix may be formed.
%       W(i).CS.rates : Sparse matrix of well rates.
%
%       W(i).CS.C     :   | MS discretization matrices wrt activeFaces
%       W(i).CS.D     :  /
%
%       W(i).CS.RHS.f :   | Coarse well system right hand sides
%       W(i).CS.RHS.h :  /
%
% SEE ALSO:
%   `computeMimeticIP`, `generateCoarseSystem`.

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


opt = struct('src', [], 'OverlapWell', 0, 'OverlapBlock', 0);
opt = merge_options(opt, varargin{:});

cellNo = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
C      = sparse  (1 : numel(cellNo), cellNo, 1);
D      = sparse  (1 : numel(cellNo), double(G.cells.faces(:,1)), 1, ...
                  numel(cellNo), G.faces.num);

% Basis weighting should always be the same as in coarse system, CS.
weight   = get_weight(G, CS, rock);
part     = get_part(CG);
numWells = numel(W);

for k = 1 : numWells,
   [Psi, Phi, R]  = evalWellBasis(W(k), G, CG, S.BI, C, D, weight, mob, ...
                                  'OverlapBlock', opt.OverlapBlock,     ...
                                  'OverlapWell' , opt.OverlapWell,      ...
                                  'src'         , opt.src);
   W(k).CS.basis  = Psi;
   W(k).CS.basisP = Phi;
   W(k).CS.rates  = R  ;

   % Determine unique coarse blocks through which the well 'W(k)' passes.
   pwc = unique(part(W(k).cells));
   ncc = numel(pwc);

   if strcmp(W(k).type, 'rate'),
      RHS.f = zeros([ncc, 1]);
      RHS.h = -W(k).val;
   else
      RHS.f = -W(k).val(ones([ncc, 1]));
      RHS.h = 0;      % never used
   end

   W(k).CS.C   = sparse(1 : ncc, pwc, 1, ncc, CG.cells.num);
   W(k).CS.D   = sparse(1 : ncc,  1 , 1, ncc,      1      );
   W(k).CS.RHS = RHS;
end

%-----------------------------------------------------------------------
% Private helpers follow
%-----------------------------------------------------------------------

function weight = get_weight(G, CS, rock)
assert (isfield(CS, 'basisWeighting'));

weight = evalBasisSource(G, CS.basisWeighting, rock);

%-----------------------------------------------------------------------

function p = get_part(cg)
p = cg.partition;
