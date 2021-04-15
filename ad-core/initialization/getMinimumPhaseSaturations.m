function s_min = getMinimumPhaseSaturations(model)
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

    nph = sum(model.getActivePhases());
    s_min = zeros(1, nph);
    
    f = model.fluid;
    sw_min = 0;
    if model.water && isfield(f, 'sWcon')
        sw_min = f.sWcon;
        ix = model.getPhaseIndex('W');
        s_min(ix) = sw_min;
    end
    
    if model.oil
        ix = model.getPhaseIndex('O');
        s_min(ix) = 0;
    end
    
    if model.gas
        ix = model.getPhaseIndex('G');
        s_min(ix) = 1 - sw_min;
    end
    
%     s = 0:0.001:1;
%     if model.water
%         krw = model.fluid.krW(s);
%         
%         ix = model.getPhaseIndex('W');        
%         s_min(ix) = getMinSat(s, krw);
%     end
%     
%     if model.oil
%         if isfield(model.fluid, 'krO')
%             kro = model.fluid.krO(s);
%         else
%             krow = model.fluid.krOW(s);
%             krog = model.fluid.krOG(s);
%             kro = min(krow, krog);
%         end
%         
%         ix = model.getPhaseIndex('O');
%         s_min(ix) = getMinSat(s, kro);
%     end
%     
%     if model.gas
%         krg = model.fluid.krG(s);
%         
%         ix = model.getPhaseIndex('G');        
%         s_min(ix) = getMinSat(s, krg);
%     end
end

function s_min = getMinSat(s, kr)
    sub = find(kr == 0, 1, 'last');
    s_min = s(sub);
end
