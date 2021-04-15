function paper = paper_ms_msrsb_polymer()
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
paper = createPaperStruct('ms-msrsb_polymer', ...
                          'Multiscale simulation of polymer flooding with shear effects', ...
                          'authors', 'S. T. Hilden, O. Moyner, K.-A. Lie and K. Bao', ...
                          'published', 'Transport in Porous Media, vol 113, issue 1, p. 111-135', ...
                          'year', 2016, ...
                          'doi', '10.1007/s11242-016-0682-2', ...
                          'modules', {'ad-eor','msrsb'}, ...
                          'url', 'http://link.springer.com/article/10.1007%2Fs11242-016-0682-2', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/mspoly.pdf');
end