function [X, Y, Z] = buildCornerPtPillars(grdecl, varargin)
%Construct physical nodal coordinates for CP grid.
%
% SYNOPSIS:
%   [X, Y, Z] = buildCornerPtPillars(grdecl)
%   [X, Y, Z] = buildCornerPtPillars(grdecl, 'pn1', pv1, ...)
%
% PARAMETERS:
%   grdecl  - Eclipse file output structure as defined by readGRDECL.
%             Must contain at least the fields 'cartDims', 'COORD' and
%             'ZCORN'.
%
% KEYWORD ARGUMENTS:
%
%  'Verbose'              - Whether or not to emit informational messages.
%                           Default value: Verbose = mrstVerbose.
%
%  'CoincidenceTolerance' - Absolute tolerance used to detect collapsed
%                           pillars where the top pillar point coincides
%                           with the bottom pillar point.  Such pillars are
%                           treated as is they were vertical.
%                           Default value: CoincidenceTolerance = `100*eps`.
%
%  'Scale'                - Scale the pillars so that the extend from the
%                           top to the bottom of the model.
%                           Default: FALSE 
%
% RETURNS:
%   X, Y, Z - Matrices with (nx+1)*(ny+1) rows with the start and end point
%             of the pillars in 'x', 'y', and 'z' direction.
%
% EXAMPLE:
%   gridfile  = [DATADIR, filesep, 'case.grdecl'];
%   grdecl    = readGRDECL(gridfile);
%   [X, Y, Z] = buildCornerPtPillars(grdecl);
%
% SEE ALSO:
%   `readGRDECL`, `processGRDECL`

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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


   opt = struct('Verbose', mrstVerbose, ...
      'CoincidenceTolerance', 100*eps, 'Scale', false);
   opt = merge_options(opt, varargin{:});

   lines = reshape(grdecl.COORD, 6, []) .';

   % Extract coordinates of pillars
   X = lines(:,[1 4]);
   Y = lines(:,[2 5]);
   Z = lines(:,[3,6]);

   % Set collapsed pillars to extend from to bottom of model
   ix = abs(lines(:,6) - lines(:,3)) < abs(opt.CoincidenceTolerance);
   zm = min(grdecl.ZCORN(:));
   zM = max(grdecl.ZCORN(:));

   Z(ix,:) = repmat([zm, zM], [sum(ix), 1]);

   if ~opt.Scale,  return,  end

   % Scale pillars to extend from top to bottom
   X(:,1) = X(:,1) + (zm - Z(:,1)).*(X(:,2)-X(:,1))./(Z(:,2)-Z(:,1));
   X(:,2) = X(:,1) + (zM - Z(:,1)).*(X(:,2)-X(:,1))./(Z(:,2)-Z(:,1));
   Y(:,1) = Y(:,1) + (zm - Z(:,1)).*(Y(:,2)-Y(:,1))./(Z(:,2)-Z(:,1));
   Z = repmat([zm zM], size(Z,1), 1);
end
