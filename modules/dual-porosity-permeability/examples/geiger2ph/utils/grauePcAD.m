function fn = grauePcAD(sr, sr_tot, model)
% Two-phase capillary pressure function handle for AD style fluids 
%
% INPUT:
%    sr - residual water saturation
%    sr_tot - sum of sr for all phases present
% DESCRIPTION:
% The capillary pressure-saturation curve is approximated as a curve fit to 
% experimental data from Graue SPE56672. The fit is based on the following data:    
%
%   sw_pc	= [0.15	0.2	0.26	0.35	0.41	0.45	0.5	0.58	0.6	0.65	0.7	0.73	0.76];
%   pc_imb_sim = [60	29	19	12	7	5	4	3	2	2	0	-1	-11] * bar_psi;
%   Swe = (sw_pc-Swr)/(1-Swr-Sor); 
%   pc = polyfit(Swe, pc_imb_sim, 3);	% -16.8901098927788	30.4949345867430	-18.0265406816741	3.78772984965509
%   pc_val = polyval(pc, Swe);

    fn = @(s) grauePc(s, sr, sr_tot, model);
end

function pc = grauePc(s, sr, sr_tot, model)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    if isa(sat, 'ADI')
        sat.val = max(min(sat.val, 1), 0);
    else
        sat = max(min(sat, 1), 0);
    end
    
    p = [-16.8901098927788, 30.4949345867430, -18.0265406816741, 3.78772984965509];   % in bars
    pc = (p(1)*sat.^3 + p(2)*sat.^2 + p(3)*sat + p(4)) * 1e5;   % in Pa
    
    % Set Pc=0 in the cells, adjacent to the outflow boundary
    %pc(model.G.iout) = 0;
    
end

%{
Copyright 2022 Geological Survey of Denmark and Greenland (GEUS).

Author: Nikolai Andrianov, nia@geus.dk.

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