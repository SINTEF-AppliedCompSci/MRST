function paper = paper_agglom_nuc()
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
paper = createPaperStruct('agglom-nuc', ...
    'Coarsening of three-dimensional structured and unstructured grids for subsurface flow', ...
    'authors', 'J. E. Aarnes, V. L. Hauge, and Y. Efendiev', ...
    'published', 'Advances in Water Resources, Vol. 30, Issue 11, pp. 2177-2193', ...
    'year', 2007, ...
    'modules', {'agglom'}, ...
    'fileurl', 'http://www.sintef.no/project/GeoScale/papers/AarnesHaugeEfendiev06.pdf', ...
    'doi', '10.1016/j.advwatres.2007.04.007', ...
    'url', 'http://www.sciencedirect.com/science/article/pii/S030917080700070X');
end
