function paper = paper_co2_fieldCase()
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
paper = createPaperStruct('co2lab-fieldCase', ...
                          'Field-case simulation of CO2-plume migration using vertical-equilibrium models', ...
                          'authors', 'H. M. Nilsen, P. A. Herrera, M. Ashraf, I. S. Ligaarden, M. Iding, C. Hermanrud, K.-A. Lie, J. M. Nordbotten, H. K. Dahle, E. Keilegavlen', ...
                          'published', 'Energy Procedia, Vol. 4, pp. 3801-3808', ...
                          'year', 2011, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/paper-ghgt.pdf', ...
                          'doi', '10.1016/j.egypro.2011.02.315', ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S1876610211005947');
end