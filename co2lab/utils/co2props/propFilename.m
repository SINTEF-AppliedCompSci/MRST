function name = propFilename(P_range, T_range, P_num, T_num, fluid, prop)
% Standardized filename generator for a sampled table of a given property of
% a given fluid, with given sample parameters.
   name = [fluid, ...
           '_', ...
           num2str(P_range, '%d_'), ...
           num2str(T_range, '%d_'), ...
           num2str(P_num,   '%d_'), ...
           num2str(T_num,   '%d_'), ...
           prop, '.mat'];
end

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