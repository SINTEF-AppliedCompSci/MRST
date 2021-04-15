function paper = paper_ms_msrsb_bo()
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
paper = createPaperStruct('ms-msrsb_bo', ...
                          'A multiscale restriction-smoothed basis method for compressible black-oil models', ...
                          'authors', 'O. Moyner and K.-A. Lie', ...
                          'published', 'SPE Journal', ...
                          'year', 2016, ...
                          'doi', '10.2118/173265-PA', ...
                          'modules', {'msfvm', 'msrsb', 'blackoil-sequential'}, ...
                          'url', 'https://www.onepetro.org/journal-paper/SPE-173265-PA', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/msrsb-bo.pdf');
end