function paper = paper_ms_fmsrsb()
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
paper = createPaperStruct('ms-fmsrsb', ...
                          'The multiscale restriction smoothed basis method for fractured porous media (F-MsRSB)', ...
                          'authors', 'S. Shah, O. Moyner, M. Tene, K.-A. Lie, and H. Hajibeygi', ...
                          'published', 'Journal of Computational Physics, Vol. 318, pp. 36-57', ...
                          'year', 2016, ...
                          'doi', '10.1016/j.jcp.2016.05.001', ...
                          'modules', {'hfm','msrsb'}, ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S0021999116301267', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/FMSRSB.pdf');
end