function [krwAr] = arithmeticAverageMPFA(G, krw, bc, pot)
% Computes the arithmetic average of the relative permeability
%
% SYNOPSIS:
%   function [krwAr] = arithmeticAverageMPFA(G, bc, krw, pot)
%
% PARAMETERS:
%   G          - Structure, Grid structure from MRST
%   bc         - Structure, Boundary conditions structure from MRST
%   krw        - Function, relative permeability function krw = krw(psi)
%                or krw = krw(p), according to the passed potential.
%   pot        - Vector, containing the values of the potential. This could
%                be either the pressure or the pressure head.
%
% RETURNS:
%   krwAr      - Vector, containing the arithmetic averaged relative
%                permeabilities at the faces.
%
% See also upstreamWeightingMPFA.

%{
Copyright 2018-2020, University of Bergen.

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

% Extracting topological data
fNei = G.faces.neighbors;               % extracting faces neighbors
int_fNei = fNei(all(fNei ~= 0,2),:);    % internal faces neighbors
int_f = find(all(fNei ~= 0,2));         % internal faces
neuFcsIdx = bc.face(all(bc.face.*strcmp(bc.type','flux') ~= 0,2),:);
dirFcsIdx = bc.face(all(bc.face.*strcmp(bc.type','pressure') ~= 0,2),:);

% Initializing output vector
krwAr = zeros(G.faces.num,1);

% Neumann boundaries
krwAr(neuFcsIdx) = 1; % equal to one, since fluxes are known

% Dirichlet boundaries
potFcsIdx = find(ismember(bc.face, dirFcsIdx));  % extracting idx
diriNeigh = fNei(dirFcsIdx,:)'; % neighboring cells of diri faces
diriCells = diriNeigh(diriNeigh > 0); % cells corresponding to each diri face
krwAr(dirFcsIdx) = 0.5.*(krw(bc.value(potFcsIdx)) + krw(pot(diriCells)));

% Internal faces
krwAr(int_f) = 0.5.*(krw(pot(int_fNei(:,1))) + krw(pot(int_fNei(:,2))));

end
