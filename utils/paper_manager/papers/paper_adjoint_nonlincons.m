function paper = paper_adjoint_nonlincons()
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
paper = createPaperStruct('adjoint_nonlincons', ...
                          'Nonlinear output constraints handling for production optimization of oil reservoirs', ...
                          'authors', 'E. Suwartadi, S. Krogstad, B. Foss', ...
                          'published', 'Computational Geosciences, Vol. 16(2), pp. 499-517', ...
                          'year', 2012, ...
                          'doi', '10.1007/s10596-011-9253-3', ...
                          'modules', {'adjoint'}, ...
                          'url', 'http://dx.doi.org/10.1007/s10596-011-9253-3');
end