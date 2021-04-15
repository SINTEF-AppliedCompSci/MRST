function paper = paper_msmfe_cmplx()
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
paper = createPaperStruct('msmfe-cmplx', ...
    'A multiscale method for modeling flow in stratigraphically complex reservoirs', ...
    'authors', 'F. O. Alpak, M. Pal, and K.-A. Lie', ...
    'published', 'SPE J., Vol. 17, No. 4, pp. 1056-1070', ...
    'year', 2012, ...
    'modules', {'msmfem'}, ...
    'doi', '10.2118/140403-PA', ...
    'url', 'https://www.onepetro.org/journals?pageType=Preview&jid=ESJ&mid=SPE-140403-PA');
end
