function paper = paper_ad_blackoilminstate()
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
paper = createPaperStruct('ad_blackoilminstate', ...
                          'Black-oil minimal fluid state parametrization for constrained reservoir control optimization', ...
                          'authors', 'A. Codas, B. Foss, E. Camponagara, S. Krogstad', ...
                          'published', 'Journal of Petroleum Science and Engineering, Vol. 143, pp. 35-43', ...
                          'year', 2016, ...
                          'doi', '10.1016/j.petrol.2016.01.034', ...
                          'modules', {'ad-core', 'ad-blackoil'}, ...
                          'url', 'http://dx.doi.org/10.1016/j.petrol.2016.01.034', ...
                          'fileurl', 'https://www.researchgate.net/profile/Andres_Codas/publication/296636014_Black-oil_minimal_fluid_state_parametrization_for_constrained_reservoir_control_optimization/links/56db08b308aebabdb412df1c.pdf');
end
