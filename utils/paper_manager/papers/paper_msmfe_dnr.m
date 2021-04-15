function paper = paper_msmfe_dnr()
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
paper = createPaperStruct('msmfe-dnr', ...
    'Grid adaptation for the Dirichlet-Neumann representation method and the multiscale mixed finite-element method', ...
    'authors', 'K.-A. Lie, J. R. Natvig, S. Krogstad, Y. Yang, and X.-H. Wu', ...
    'published', ' Comput. Geosci., Vol. 18, No. 3, pp. 357-372', ...
    'year', 2014, ...
    'modules', {'msmfem'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/dnr-comg.pdf', ...
    'doi', '10.1007/s10596-013-9397-4', ...
    'url', 'http://link.springer.com/article/10.1007%2Fs10596-013-9397-4');
end
