function refSol = eclipseOutput2Sol(G, output, active)
%Convert ECLIPSE output to MRST rSol object.
%
% SYNOPSIS:
%   refSol = eclipseOutput2Sol(G, output, active)
%
% PARAMETERS:
%   G      - Grid data structure.  Must contain valid field
%            'G.cells.indexMap'.
%
%   output - ECLIPSE output structure as defined by function
%            'readEclipseOutput'.
%
%   active - List of explicitly active cells.  These are the global (i.e.,
%            linearised Cartesian) cells for which ACTNUM==1 and, if
%            applicable, PORO>0.
%            Logical array of size PROD(G.cartDims)-by-1, the TRUE elements
%            of which correspond to the explicitly active cells.
%
% RETURNS:
%   refSol - Reference solution as an MRST rSol-type object/structure.

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


   [sat, pres] = deal(repmat(-1, [prod(G.cartDims), 1]));

   assert (numel(output.SWAT.values) == numel(output.PRESSURE.values));
   assert (numel(output.SWAT.values) == sum(double(active)));

   sat (active) = output.SWAT.values;
   pres(active) = output.PRESSURE.values;
   refSol.s     = sat(G.cells.indexMap);
   refSol.pressure = convertFrom(pres(G.cells.indexMap), barsa);
end
