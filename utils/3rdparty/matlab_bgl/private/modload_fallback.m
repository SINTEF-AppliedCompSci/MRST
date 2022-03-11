function modload
%Fallback Strategy for MATLAB BGL Activation Failure

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

   d = fileparts(fileparts(mfilename('fullpath')));

   error(['MatlabBGL software is not available\n\n',    ...
          ' -> Use downloadMBGL script to install ',    ...
          '(copy-paste command below)\n\n',             ...
          '    run ''%s''\n\n',                         ...
          'Then run mrstModule add matlab_bgl again.'], ...
          fullfile(d, 'downloadMBGL'));
end
