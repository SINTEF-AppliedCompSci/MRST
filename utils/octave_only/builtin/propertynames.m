function v = propertynames(x)
%Undocumented Utility Function

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

  % disable warning in call to struct(x) below
  warn_id = 'Octave:classdef-to-struct';
  initial_warn_id_state = warning('query', warn_id).state;

  warning('off', warn_id);
  
  all_props = fieldnames(struct(x)); % includes public, protected and
                                     % private properties

  warning(initial_warn_id_state, warn_id);
  
  % keep only public properties
  v = {};
  for i = 1:numel(all_props)
     name = all_props{i};
     try
        tmp = x.(name);
      
        % this method was accessible, so keep it
        v = [v; {name}];
     catch
      
        % this method was inaccessible, i.e. not public.  Discard it
     end_try_catch
  end
end
