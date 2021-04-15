function paper = paper_co2_calibrateSleipner()
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
paper = createPaperStruct('co2lab-calibrateSleipner', ...
                          'Using Sensitivities and Vertical-equilibrium Models for Parameter Estimation of CO2 Injection Models with Application to Sleipner Data', ...
                          'authors', 'H. M. Nilsen, S. Krogstad, O. Andersen, R. Allen, K.-A. Lie', ...
                          'published', 'Energy Procedia, Vol. 114, pp. 3476-3495', ...
                          'year', 2017, ...
                          'modules', {'co2lab'}, ...
                          'fileurl','http://folk.ntnu.no/andreas/papers/Sleipner_GHGT.pdf', ...
                          'doi', '10.1016/j.egypro.2017.03.1478', ...
                          'url', 'https://www.sciencedirect.com/science/article/pii/S1876610217316727');
end