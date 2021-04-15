function paper = paper_co2_spillPoint()
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
paper = createPaperStruct('co2lab-spillPoint', ...
                          'Spill-point analysis and structural trapping capacity in saline aquifers using MRST-co2lab', ...
                          'authors', 'H. M. Nilsen, K.-A. Lie, O. Moyner, O. Andersen', ...
                          'published', 'Computers and Geoscience, Vol. 75, pp. 33-43', ...
                          'year', 2015, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/co2lab-1.pdf', ...
                          'doi', '10.1016/j.cageo.2014.11.002', ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S0098300414002520');
end