function [pcOW, pcOG] = getCapillaryPressureBO(f, sW, sG)
    [pcOW, pcOG] = deal(0);
    
    % Check for capillary pressure (p_cow)
    if isfield(f, 'pcOW') && ~isempty(sW)
        pcOW  = f.pcOW(sW);
    end
    
    % Check for capillary pressure (p_cog)
    if isfield(f, 'pcOG') && ~isempty(sG) && nargin > 1
        pcOG  = f.pcOG(sG);
    end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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
