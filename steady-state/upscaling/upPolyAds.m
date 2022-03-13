function updata = upPolyAds(block, updata, varargin)
% Upscale polymer adsorption isotherm using a simple average.
% 
% The adsorption isotherm which is a function of the local polymer solution
% concentration.

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

opt = struct(...
    'npoly', 50 ...
    );
opt = merge_options(opt, varargin{:});

G = block.G;
f = block.fluid;
c = linspace(0, f.cmax, opt.npoly)';

% The adsorption is computed by a rock mass average
rm    = f.rhoR .* G.cells.volumes .* (1 - block.rock.poro);
rmtot = sum(rm);
ads   = nan(opt.npoly, 1);
for i=1:opt.npoly
   ads(i) = sum(rm .* f.ads( c(i).*ones(G.cells.num,1) ) ) / rmtot;
end

% Add to updata structure
updata.ads = [c ads];

end
