struct_levels_to_print(0); % nicer output of nested structs

% remove recurring warning triggered by out workaround to Octaves
% incomplete support for classes
warn_id = 'Octave:classdef-to-struct';
warning('off', warn_id);
warning('off', 'Octave:lu:sparse_input')

page_output_immediately (1);
page_screen_output (0);
warning('off', 'Octave:divide-by-zero')

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
