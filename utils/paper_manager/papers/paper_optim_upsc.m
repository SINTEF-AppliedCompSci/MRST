function paper = paper_optim_upsc()
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
paper = createPaperStruct('optim_upsc', ...
                          'Reservoir management optimization using well-specific upscaling and control switching', ...
                          'authors', 'S. Krogstad, X. Raynaud, H.M. Nilsen', ...
                          'published', 'Computational Geosciences, Vol. 20(3), pp. 695-706', ...
                          'year', 2016, ...
                          'doi', '10.1007/s10596-015-9497-4', ...
                          'modules', {'optimization','upscaling'}, ...
                          'url', 'http://dx.doi.org/10.1007/s10596-015-9497-4');
end