function [pc_upsc, pc_max, pc_min, sat_min, sat_max] = upscalePC(G, fluid, pv)
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


    pc_min=min(fluid.pc(ones(G.cells.num,1)),1);
    pc_max=man(fluid.pc(zeros(G.cells.num,1)),1);
    pc_val=linespace(pc_min,pc_max,100);
    volume=sum(G.cells.volumes.*pv);
    sat_min=1;
    sat_max=0;
    for i=1:numel(pc_val)
        s_cells=fluid.invpc(pc_val*ones(G.cells.num,1));
        s_val=sum(s_cells.*pv)/volume;
    end
    pc_upsc=@(s) pcUpscaled(s, s_val, pc_val)

end
function varargout = pcUpcaled(s, s_val, pc_val)
  varargout{1}    = interpTable(s_val, pc_val, s);

  if nargout > 1,
     varargout{2} = dinterpTable(s_val, pc_val, sg);
  end
end
