function min_rs= minRs(p,sG,sGmax,f, G)
% Determine the minimal amount of dissolved CO2 in each cell, based on the
% maximum historical CO2 saturation in each cell.  The brine in the part of
% the column that contains CO2 (whether free-flowing or residual) is assumed
% to be saturated with CO2.
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
    % computing |g|.(rho_w - rho_g)
    drho=norm(gravity)*(f.rhoWS.*f.bW(p)-f.rhoGS.*f.bG(p));
    
    % Computing position of the h_max interface, consistent with the current
    % value of sGmax.
    pcmax=f.pcWG(sGmax, p,'sGmax',value(sGmax));       
    h_max=pcmax./drho;
    assert(all(value(h_max)>=0));

    % If the computed value for h_max ever exceeds H, set it to H (while
    % keeping ADI structure if applicable).
    ind=value(h_max)>G.cells.H;
    h_max(ind)=G.cells.H(ind)+h_max(ind)*0;
    
    % Minimal dissolved quantity equals the 'dis_max' value multiplied by the
    % content of brine in the zone above 'h_max' (i.e. residual brine within the
    % CO2 plume, as well as the brine in the residual CO2 zone).  The content of
    % brine in this zone equals the total content of brine in the column (1-sG),
    % minus the content of brine below 'h_max'.
    min_rs=((1-sG)-(G.cells.H-h_max)./G.cells.H).*f.dis_max;   

end

