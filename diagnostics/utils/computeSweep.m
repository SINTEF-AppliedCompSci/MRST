function [Ev,tD] = computeSweep(F,Phi)
%Compute sweep efficiency versus dimensionless time (PVI)
%
% SYNOPSIS:
%   [Ev,tD] = computeSweep(F, Phi)
%
% PARAMETERS:
%   F   - flow capacity
%   Phi - storage capacity
%
% RETURNS:
%   tD  - dimensionless time (PVI)
%   Ev  - sweep efficiency. Analogue: 1D displacement using the F-Phi curve
%         as flux funcion
%
% SEE ALSO:
%   `computeFandPhi`, `computeLorenz`

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


% remove any null segments (to avoid division by zero)
hit        = true(size(F));
hit(2:end) = F(1:end-1) <= F(2:end) - sqrt(eps);
F          = F(hit);
Phi        = Phi(hit);

% compute dimensionless time
tD = [0; (Phi(2:end)-Phi(1:end-1))./(F(2:end)-F(1:end-1))];

% compute sweep
Ev   = Phi + (1-F).*tD;

end

