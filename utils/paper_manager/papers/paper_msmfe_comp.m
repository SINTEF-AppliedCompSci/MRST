function paper = paper_msmfe_comp()
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
paper = createPaperStruct('msmfme-comp', ...
    'A comparison of multiscale methods for elliptic problems in porous media flow', ...
    'authors', 'V. Kippe, J. E. Aarnes, and K.-A. Lie', ...
    'published', 'Comput. Geosci. (Special issue on multiscale methods), Vol. 12, No. 3, pp. 377-398', ...
    'year', 2008, ...
    'modules', {'msmfem'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/ms-compare.pdf', ...
    'doi', '10.1007/s10596-007-9074-6', ...
    'url', 'http://link.springer.com/article/10.1007%2Fs10596-007-9074-6');
end
