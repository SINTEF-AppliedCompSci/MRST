function [pc_upsc, pc_max, pc_min, sat_min, sat_max] = upscalePCCaplimit(G, fluid, pv)
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


    % hack to find maximum/minimum capillary pressure
    % assumes finescale pc gives the maximum/minimum out of range
    state=struct('s',ones(G.cells.num,1));
    state.extSat=[state.s,1-state.s];
    pc_min=min(fluid.pc(state));
    state.s=zeros(G.cells.num,1);
    state.extSat=[state.s,1-state.s];
    pc_max=max(fluid.pc(state));

    if(pc_min<pc_max)
        pc_val=linspace(pc_min,pc_max,100);
        volume=sum(G.cells.volumes.*pv);
        %sat_min=1;
        %sat_max=0;
        s_val=nan(numel(pc_val),1);
        for i=1:numel(pc_val)
           s_cells=fluid.invpc(pc_val(i)*ones(G.cells.num,1));
          s_val(i)=sum(s_cells.*G.cells.volumes.*pv)/volume;
        end
        sat_min=s_val(1);
        sat_max=s_val(end);
        pc_upsc=@(s) pcUpscaled(s, s_val, pc_val);
    else
       sat_min=1;sat_max=0;
       pc_upsc=@(s) ones(size(s))*pc_min;
    end
end
function varargout = pcUpscaled(s, s_val, pc_val)
  varargout{1}    = interpTable(s_val, pc_val, s);

  if nargout > 1,
     varargout{2} = dinterpTable(s_val, pc_val, sg);
  end
end
