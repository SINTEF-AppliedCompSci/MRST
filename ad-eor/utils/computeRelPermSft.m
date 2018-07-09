function [krW, krO] = computeRelPermSft(sW, c, Nc, fluid)
%
%
% SYNOPSIS:
%   function [krW, krO] = computeRelPermSft(sW, c, Nc, fluid)
%
% DESCRIPTION: Computes the water and oil relative permeabilities, using the
% surfactant model as described in  ad-eor/docs/surtactant_model.pdf
%
% PARAMETERS:
%   sW    - Saturation
%   c     - Concentration
%   Nc    - Capillary number
%   fluid - Fluid structure
%
% RETURNS:
%   krW - Water relative permeability
%   krO - Oil relative permeability
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright 2009-2018 SINTEF ICT, Applied Mathematics.

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


    isSft = (double(c) > 0);
    m = 0*c;
    if nnz(isSft) > 0
       logNc = log(Nc(isSft))/log(10);
       % We cap logNc (as done in Eclipse)
       logNc = min(max(-20, logNc), 20);
       m(isSft) = fluid.miscfact(logNc, 'cellInx', find(isSft));
    end

    sWcon    = fluid.sWcon;    % Residual water saturation   without surfactant
    sOres    = fluid.sOres;    % Residual oil saturation     without surfactant
    sWconSft = fluid.sWconSft; % Residual water saturation   with    surfactant
    sOresSft = fluid.sOresSft; % Residual oil saturation     with    surfactant

    % Interpolated water/oil residual saturations
    sNcWcon = m.*sWconSft + (1 - m).*sWcon;
    sNcOres = m.*sOresSft + (1 - m).*sOres;

    sNcEff = (sW - sNcWcon)./(1 - sNcWcon - sNcOres);

    % Rescaling of the saturation - without surfactant
    sNcWnoSft = (1 - sWcon - sOres).*sNcEff + sWcon;
    krNcWnoSft = fluid.krW(sNcWnoSft);
    krNcOnoSft = fluid.krOW(1 - sNcWnoSft);

    % Rescaling of the saturation - with surfactant
    sNcWSft =  (1 - sWconSft - sOresSft).*sNcEff + sWconSft;
    krNcWSft = fluid.krW(sNcWSft);
    krNcOSft = fluid.krOW(1 - sNcWSft);

    krW = m.*krNcWSft + (1 - m).*krNcWnoSft;
    krO = m.*krNcOSft + (1 - m).*krNcOnoSft;


end
