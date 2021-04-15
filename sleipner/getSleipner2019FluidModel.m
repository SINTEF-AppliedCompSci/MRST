function [fluid] = getSleipner2019FluidModel(varargin)
% Function to get fluid properties for the Sleipner 2019 benchmark
%
% SYNOPSIS:
%   [fluid,tsurf,tgrad,pref,rhoCref] = getSleipner2019FluidModel(varargin)
%
% OPTIONAL PARAMETERS:
%
%   lincomp     - Flag to turn on linear compressibility for the gas phase.
%               This is only relevant when not using the EOS.
%   n           - relperm exponents
%   useEOS      - Flag to use coolprop EOS to determine CO2 density 
%   topReservoirTemp - Temperature at top of the reservoir.
%   injTemp     - Temperature at injection point
%   G           - Grid structure, required when using EOS. 
%   wcell       - Cell containing injection well. Used to determine depth 
%                    of injection point. Default is injection cell in 
% 
% RETURNS:
%   fluid - fluid structure with Sleipner specific parameters.

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


opt = struct('lincomp', true,...       % use constant compressibility co2
            'n',[1 1],...
             'useEOS',true, ...
             'topReservoirTemp', 37, ...
             'injTemp',48, ...
             'G',[], ...
             'wcell',70645); % This is the first well cell in the multilayeredVE model

opt = merge_options(opt, varargin{:});

gravity reset on

%% Specify fluid properties:
% Here, we use the benchmark values except for the residual
% saturations. Seafloor depth info is not given in benchmark, thus a
% depth of 100m is assumed. Compressibilities were calculated using
% Span and Wagner's EOS (via coolprops) for reservoir conditions
% corresponding to initial pressure and temperature (i.e., hydrostatic
% and computed via thermal gradient), however we set CO2
% compressibility to be 10 times that of water.
[rho, mu]   = getValuesSPE134891();

water_density       = rho(1) * kilogram / meter ^3;
rhoCref             = rho(2) * kilogram / meter ^3;


water_compr_val     = 4.37e-5/barsa; % will convert to compr/Pa
if(opt.lincomp)
    co2_compr_val       = water_compr_val*10; % 1.24e-3/barsa;
else
    co2_compr_val   = [];
end

pref               = 6*mega*Pascal;             % ref. for linear compressibilities


%% Create fluid:
fluid = initSimpleADIFluid('phases', 'WG'           , ...
    'mu'  , [mu(1), mu(2)]     , ...
    'rho' , [water_density, rhoCref]     , ...
    'pRef', pref           , ...
    'c'   , [water_compr_val co2_compr_val],...
    'n'   , opt.n);



if opt.useEOS
    injT = opt.injTemp; % Celsius
    topT = opt.topReservoirTemp; % Celsius
    assert(~isempty(opt.G),'Need G to calc. temp. dependent density');
    G = opt.G;
    injDepth = G.cells.centroids(opt.wcell,3);
    topDepth = min(G.cells.centroids(:,3));
    gradT = (injT - topT)./(injDepth - topDepth);
    
    reservoir_temp = 273.15 + topT + ...
        (G.cells.centroids(:,3) - topDepth) * gradT;    % Kelvin

    fluid = addSampledFluidProperties(fluid,...
                                      'G','fixedT',reservoir_temp);

    fluid = include_BO_form(fluid, 'G', fluid.rhoGS);
end

    

