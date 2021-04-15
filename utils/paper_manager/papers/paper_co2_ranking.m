function paper = paper_co2_ranking()
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
paper = createPaperStruct('co2lab-ranking', ...
                          'Ranking and categorizing large-scale saline aquifer formations based on optimized CO2 storage potentials and economic factors', ...
                          'authors', 'R. Allen, H. M. Nilsen, O. Andersen, K.-A. Lie', ...
                          'published', 'International Journal of Greenhouse Gas Control, Vol. 65, pp. 182-194', ...
                          'year', 2017, ...
                          'modules', {'co2lab'}, ...
                          'fileurl','https://folk.ntnu.no/andreas/papers/IJGGC-paper.pdf', ...
                          'doi', '10.1016/j.ijggc.2017.07.023', ...
                          'url', 'https://www.sciencedirect.com/science/article/pii/S1750583617302657?via%3Dihub');
end