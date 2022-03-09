function fluid = addThermalFluidProps( fluid, varargin )
%Add thermal properties to an existing fluid structure
%
% SYNOPSIS:
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, ...);
%  fluid = addThermalFluidProps(fluid,'pn1', pv1, 'useEOS', true, 'brine', true, ...);
%
% PARAMETERS:
%   fluid   - fluid structure created with initSimpleADIFluid (or other).
%
%   cp      - fluid heat capacity in J kg-1 K-1. typically 4.2e3. Used for
%             the evaluation of the internal energy and enthalpy.
%
%   lambdaF - fluid heat conductivity in W m-1 K-1. typically 0.6.
%
%   useEOS  - logical. By default the density/viscosity formulation
%              of Spivey is used.
%
%   brine   - logical. Used for the EOS and p,T,c dependency of properties
%             if salt is present or not.
%
%   dNaCl   - salt molecular diffusivity in m2s-1.
% 
% 
%   rho     - density at reservoir condition. Can be given as an handle 
%             function.
% 
%   useBFactor - logical. Used if we want bW definition factor instead 
%                of rhoW 
% 
% RETURNS:
%   fluid - updated fluid structure containing the following functions
%           and properties.
%             * bX(p,T,c)  - inverse formation volume factor
%             * muX(p,T,c) - viscosity functions (constant)
%             * uW(p,T,c)  - internal energy
%             * hW(p,T,c)  - enthalpy
%             * lambdaF    - fluid conductivity
%             * dNaCl      - Salt molecular diffusivity
%
%
% SEE ALSO:
%
% 'initSimpleADIFluid', 'initSimpleThermalADIFluid'

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

    Watt = joule/second;
    opt = struct('Cp'     , 4.2*joule/(Kelvin*gram), ...
                 'lambdaF', 0.6*Watt/(meter*Kelvin), ...
                 'useEOS' , false                  , ...
                 'rho'    , fluid.rhoWS            , ...
                 'cT'     , []                     , ...
                 'TRef'   , (273.15 + 20)*Kelvin   , ...
                 'cX'     , []                     );
    opt = merge_options(opt, varargin{:});
    % Set thermal conductivity
    fluid.lambdaF = opt.lambdaF;
    % Get phase names
    fNames = fieldnames(fluid);
    phases = '';
    for fNo = 1:numel(fNames)
        fn = fNames{fNo};
        if strcmpi(fn(1:2), 'mu')
            phases = [phases, fn(3)];
        end
    end
    % Set viscosity and density
    nPh = numel(phases);
    if any(strcmpi(phases, 'W')) && opt.useEOS
        fluid.muW  = @(varargin) computeEOSViscosity(varargin{:});
        fluid.rhoW = @(varargin) computeEOSDensity(varargin{:});
        fluid.rhoWS = fluid.rhoW(1*atm, opt.TRef, 0);
    else
        fluid.rhoW = @(varargin) computeSimpleDensity(fluid, opt, varargin{:});
    end
    % Set specific heat capacity, internal energy and enthalpy
    names = upper(phases);
    Cp    = opt.Cp;
    for phNo = 1:nPh
        n = names(phNo);
        fluid.(['Cp', n]) = @(varargin) Cp(phNo);
        fluid.(['u', n])  = @(varargin) computeInternalEnergy(fluid, n, varargin{:});
        fluid.(['h', n])  = @(varargin) computeEnthalpy(fluid, n, varargin{:});
    end
end

%-------------------------------------------------------------------------%
function mu = computeEOSViscosity(varargin)
    p = varargin{1};
    T = varargin{2};
    if nargin == 2 || isempty(varargin{3})
        mu = viscosity_pure_water(p, T);
    else
        X = varargin{3};
        mu = viscosity_brine(p, T, X);
    end
end

%-------------------------------------------------------------------------%
function rho = computeEOSDensity(varargin)
    p = varargin{1};
    T = varargin{2};
    if nargin == 2 || isempty(varargin{3})
        rho = density_pure_water(p, T);
    else
        X = varargin{3};
        rho = density_brine(p, T, X);
    end
end

%-------------------------------------------------------------------------%
function rho = computeSimpleDensity(fluid, opt, varargin)
    p = varargin{1};
    T = varargin{2};
    if nargin < 6, X = 0; else, X = varargin{3}; end
    c = fluid.bW(p);
    if ~isempty(opt.cT)
        c = c.*exp(-opt.cT*(T - opt.TRef));
    end
    if ~isempty(opt.cX)
        c = c.*exp(opt.cX*X);
    end
    rho = fluid.rhoWS.*c;
end

%-------------------------------------------------------------------------%
function u = computeInternalEnergy(fluid, phase, varargin)
    Cp = feval(fluid.(['Cp', phase]), varargin{:});
    T  = varargin{2};
    u  = Cp*T;
end

%-------------------------------------------------------------------------%
function h = computeEnthalpy(fluid, phase, varargin)
    rho = feval(fluid.(['rho', phase]), varargin{:});
    u   = feval(fluid.(['u', phase]), varargin{:});
    p   = varargin{1};
    h   = u + p./rho;
end