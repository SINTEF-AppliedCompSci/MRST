function [krwUp] = upstreamWeightingMPFA(G, bc, bcVal, mpfa_discr, ...
    phys, pot, typePot, grav)
% Computes the upstream weighting of the relative permeability
%
% SYNOPSIS:
%   function krwUp = upstreamWeightingMPFA(G, bc, bcVal, mpfa_discr, krw, gamma, pot, typePot, grav)
%
% PARAMETERS:
%   G           - Structure, Grid structure from MRST.
%   bc          - Structure, Boundary conditions structure from MRST
%   bcVal       - Vector, containing the values of boundary conditions.
%                 In the case Dirichlet conditions are prescribed and
%                 gravity effects considered, bcVal is composed by the 
%                 total contribution from both.
%   mpfa_discr  - Structure, MPFA discretization structure.
%   krw         - Function, relative permeability function krw = krw(psi)
%                 or krw = krw(p), according to the passed potential.
%   gamma       - Scalar, specific gravity, i.e., gamma = rho * g
%   pot         - Vector, containing the values of the potential. 
%   typePot     - String, type of potential. This could be either 'psi'
%                 for pressure head, or 'pressure' for pressure.
%   grav        - String, 'on' when gravity effects are included in the 
%                 model, and 'off' when they are not.
%
%  RETURNS:
%   krwUp       - Vector, containing the upstream weighted relative
%                 permeabilities at the faces. 

%{
Copyright 2018-2019, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%} 

% Retrieving relative permeability from SWRC
[~, krw, ~] = vanGenuchtenMualemSw(phys.flow.a, phys.flow.S_r, ...
    phys.flow.n, phys.flow.m);

% Extracting grid and boundary information
fNei = G.faces.neighbors;               % extracting faces neighbors
int_fNei = fNei(all(fNei ~= 0,2),:);    % internal faces neighbors
int_f = find(all(fNei ~= 0,2));         % internal faces
neuFcsIdx = bc.face(all(bc.face .* strcmp(bc.type', 'flux') ~= 0,2),:);
dirFcsIdx = bc.face(all(bc.face .* strcmp(bc.type', 'pressure') ~= 0,2),:);
krwUp = zeros(G.faces.num,1);           % initializing the upstream krw

% Obtaining elevation head
zetac = max(G.faces.centroids(:, end)) - G.cells.centroids(:, end);

% Neumann boundaries relative permeabilities
krwUp(neuFcsIdx) = 1;   % Since the flux is already known, we set krw=1

% Computing the flux at every face of the grid
F = @(x) mpfa_discr.F * x;
boundF = @(x) mpfa_discr.boundFlux * x;

% Check if gravity should be included
if strcmp(grav, 'off')
    zetac(:) = 0;
elseif ~strcmp(grav, 'on')
    error('Gravity argument not recognized. Use either ''on'' or ''off''')
end

% Compute fluxes checking if psihead or pressure is used. We assume that
% bcVal contains the gravity contribution, in the case grav = 'on'.
if strcmp(typePot, 'psiHead')
    flux = F(pot + zetac) + boundF(bcVal);
elseif strcmp(typePot, 'pressure')
    flux = F(pot + phys.flow.gamma.*zetac) + boundF(bcVal);
else
    error('Potential not recognized. See documentation for valid potentials.')
end

% Dirichlet boundaries relative permeabilities
potFcsIdx = find(ismember(bc.face,dirFcsIdx));  % extracting h at every dirichlet face from the bc struct
diriNeigh = fNei(dirFcsIdx,:)'; % neighboring cells of dirichlet faces
diriCells = diriNeigh(diriNeigh > 0); % cells corresponding to each dirichlet face
upIsFace = fNei(dirFcsIdx, 1) == 0;   % is the upstream direction a face? ~upIsFace corresponds to a cell
downIsFace = fNei(dirFcsIdx, 2) == 0; % is the downstream direction a face? ~downIsFace corresponds to a cell
% We compute the flux at each dirichlet face. If the flux >= 0 then we
% evaluate in the upstream direction. Then, we ask if the upstream
% direction is a face or not. If true, we evaluate at the face, if not
% at the cell center. If the flux < 0, we evaluate in the downstream
% direction. Then, we ask if the downstream direction is a face or not.
% If true,  we evaluate at the face, if not at the cell center.

krwUp(dirFcsIdx) = ( ...
    (flux(dirFcsIdx) > 0) .* upIsFace    .* krw(bc.value(potFcsIdx)) + ...
    (flux(dirFcsIdx) > 0) .* ~upIsFace   .* krw(pot(diriCells))      + ...
    (flux(dirFcsIdx) <= 0)  .* downIsFace  .* krw(bc.value(potFcsIdx)) + ...
    (flux(dirFcsIdx) <= 0)  .* ~downIsFace .* krw(pot(diriCells)) ...
    );

% Internal faces
% upCenters determines the indices at which the relative permeabilities must
% be evaluated. If flux >= 0 then krw is evaluated at the upstream
% direction. If flux < 0 then krw is evaluated at the downstream direction.

upCenters = ( ...
    (flux(int_f) > 0) .* int_fNei(:,1) + ... % evaluate at i-1
    (flux(int_f) <= 0)  .* int_fNei(:,2)   ... % evaluate at i
    );
krwUp(int_f) = krw(pot(upCenters));

end
