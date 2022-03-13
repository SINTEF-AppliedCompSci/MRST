function [sgn, coarse_f] = signOfFineFacesOnCoarseFaces(G, CG, nsub, sub)
%Identify fine-scale flux direction corresponding to coarse-scale outflow
%
% SYNOPSIS:
%   [sgn, coarse_f] = signOfFineFacesOnCoarseFaces(G, CG, nsub, sub)
%
% DESCRIPTION:
%   Generates a vector sgn which gives the sign that the fine fluxes should
%   be multiplied with, before summed along a coarse grid interface to
%   obtain a coarse net flux. The vector coarse_f gives the coarse face
%   number corresponding to the fine grid faces in the coarse grid.
%
% REQUIRED PARAMETERS:
%   G    - Grid structure
%
%   CG   - Coarse grid sturcture.
%
%   nsub - Number of fine-grid faces belonging to each individual coarse
%          grid face.
%
%   sub  - Actual fine-grid subfaces represented in a condensed storage
%          format.
%
% RETURNS:
%   sgn      - Sign of the fine-grid subface flux to be multiplied with
%              the corresponding flux, before summed along the coare
%              interface to give the coarse net flux.
%
%   coarse_f - The coarse face number for each fine-grid subface.
%
% SEE ALSO:
%   `subFaces`

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


   coarse_f    = rldecode(1:numel(nsub), nsub, 2) .';
   outBlckC    = rldecode(CG.faces.neighbors(:,1), nsub);

   p           = [0; CG.partition];
   outBlckF    = p(G.faces.neighbors(sub,1) + 1);

   sgn         = 2*(outBlckF == outBlckC) - 1;
end
