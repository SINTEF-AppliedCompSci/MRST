function paper = paper_co2_trapCap()
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
paper = createPaperStruct('co2lab-trapCap', ...
                          'Analysis of CO2 trapping capacities and long-term migration for geological formations in the Norwegian North Sea using MRST-co2lab', ...
                          'authors', 'H. M. Nilsen, K.-A. Lie, O. Andersen', ...
                          'published', 'Computers and Geoscience, Vol. 79, pp. 15-26', ...
                          'year', 2015, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/co2lab-4.pdf', ...
                          'doi', '10.1016/j.cageo.2015.03.001', ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S0098300415000497');
end