function x = initResSolComp(G, W, fluid, initP, massfraction)
%Initialise a state object using compressible fluids
%
% SYNOPSIS:
%   state = initResSolComp(G, W, fluid, initP, massfrac)
%
% PARAMETERS:
%   G, W     - Grid and well structures.  If the well structure 'W' is
%              empty the function will create a state object that does not
%              account for wells.
%
%   fluid    - Compressible fluid data structure.
%
%   initP    - Initial reservoir pressure.
%
%   massfrac - Mass fraction of reservoir fluid components.
%
% RETURNS:
%   state - Fully initialised reservoir/well solution structure.
%
% SEE ALSO:
%   initSingleCompFluid.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   if ischar(massfraction),
      switch massfraction,
         case 'saturated oil',
            p = initP; z0 = [1, 1, inf]';

            for i = 1:10,
               [u, u, u, u, R, B, A] = fluid.pvt(p, z0');              %#ok
               z0 = A * [0; 1; 0];
            end

            z0 = reshape(z0, 1, []);
            ig = strcmpi('gas', fluid.names);

            z0(ig) = z0(ig) * 10;

            if ~ (u(ig) < sqrt(eps)),
               warning('Ugas >= SQRT(eps)');
            end

         otherwise,
            error('''%s'' mode not supported', massfraction);
      end
   else
      assert (isnumeric(massfraction), 'Huh?');
      z0    = bsxfun(@rdivide, massfraction, fluid.surfaceDensity);
   end

   [u, u, u, u] = fluid.pvt(initP, z0);                                %#ok
   alpha = sum(u, 2);
   z0    = bsxfun(@rdivide, z0, alpha);
   s0    = bsxfun(@rdivide, u,  alpha);

   x     = initResSol (G, initP, s0, z0);

   % Set plausible initial face pressure
   cellno = rldecode(1:G.cells.num, diff(G.cells.facePos), 2).';
   x.facePressure = ...
               accumarray(G.cells.faces(:,1), ...
                          x.pressure(cellno), [G.faces.num, 1]) ./ ...
               accumarray(G.cells.faces(:,1), 1, [G.faces.num, 1]);

   if ~isempty(W),
      x.wellSol = initWellSol(W, initP);
   else
      x.wellSol = [];
   end
end
