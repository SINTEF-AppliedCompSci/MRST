function fsg = free_sg(sg, sGmax, opt)
% Determine the mobile part of present saturation.
% 
% SYNOPSIS:
%   function fsg = free_sg(sg, sGmax, opt)
%
% DESCRIPTION:
% Assuming a sharp interface, this function determine the amount of present
% saturation that is in the 'mobile plume' domain (as opposed to the region
% below the mobile plume where CO2 has been residually trapped after
% imbibition.  As such,  the "mobile part of present saturation" does include
% the CO2 in the mobile plume that is destined to be left behind as residual
% trapping, but not the CO2 that is already residually trapped.
% 
% The formula is based on the simple transformation: 
%     s * H = h * (1 - sr(2)) + (h_max - h) * sr(1)
% s_max * H = h_max * (1 - sr(2))
% 
% PARAMETERS:
%   sg    - present saturation
%   sGmax - historically maximum saturation
%   opt   - structure expected to contain the following fields:
%           * opt.res_gas : residual gas saturation
%           * opt.res_water : residual oil saturation
% RETURNS:
%   fsg - the free part of the present saturation (fsg <= sg)
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

    %% aliases for better readability
    rw = opt.res_water; % 'wetting phase'
    rn = opt.res_gas; % 'non-wetting phase'
    
    %% Computing free part of current saturation
    % NB: in the h-formulation, this would equal (1-sw) h / H. 

    fsg = ((1 - rw) * sg - (sGmax * rn)) ./ (1 - rw - rn);

    %% Ensuring exact values at boundaries
    ineb = (sg >= sGmax);  % cells with no residual trapping (drainage zone)
    fsg = ifcond(sg, fsg, ineb);
    fsg = max(fsg, 0.0);
end 
