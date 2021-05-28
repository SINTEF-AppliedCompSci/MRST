function fluid = initSimpleADIFluid(varargin)
% Make a structure representing a fluid with up to three phases (water,
% oil, gas) and their properties (relative permeabilities, densities,
% viscosities). Relative permeabilities are modeled as monomial functions
% of saturation. Densities are assumed constant, so the returned formation
% volume factors are constant equal to 1, unless compressibility factors
% are given as optional arguments.
%
% SYNOPSIS:
%   fluid = initSimpleADIFluid()
%   fluid = initSimpleADIFluid('pn1', 'pv1')
%
% OPTIONAL PARAMETERS:
%   phases - A string containing up to three of the letters W, O, G,
%            representing water, oil and gas respectively. The values of
%            this input argument determines the interpretation of the other
%            input arguments. For instance, the default is WOG, which will
%            lead to the 'rho' input argument
%            ... 'rho', [1000, 700, 100]*kilogram/meter^3)
%            to be interpreted as a water phase with density of 1000kg/m^3,
%            the oil phase having a density 700 kg/m^3 and gas having a
%            density of 100. However, if the input argument was 
%            ... 'rho', [300, 500])
%            and 'phases' was 'GO', it would be interpreted as a gas-oil
%            system with gas density of 300kg/m^3 and oil density of 500
%            kg/m^3. See the examples and the end of the documentation for
%            more information.
%
%   mu    - vector of viscosity values for phases present. For a
%           three-phase model, this will be [muW, muO, muG].
%           Default is [1 1 1].
%
%   rho   - vector of surface density values for phases present. For a
%           three-phase model, this will be [rhoWS, rhoOS, rhoGS].
%           Default is [1 1 1].
%
%   n    - vector of the degrees of the monomials describing relative
%           permeability for each phase. Default is [1 1 1] (linear
%           relative permeabilities for each phase). If the model is a
%           three-phase model, oil-water and oil-gas relative permeability
%           curves will also be added.
%
%   b     - Inverse formation volume factors. The reservoir density of a
%           phase is defined as rho_reservoir = rho_surf * b. The inverse
%           formation volume factor, denoted small b, is the reciprocal of
%           the large B often seen in literature,
%           B = 1/b.
%
%   c     - Compressibility factor. If specified for each phase, it will
%           result in b-factors on the form ::
%
%               b(p) = b_ref * exp((p-pRef)*c),
%
%           where b_ref is the b-factor specified as keyword argument (see
%           above) and pRef is specified as a seperate keyword (default: 0)
%
%   cR    - Rock compressibility. If given, the fluid will contain a rock
%           pore volume multiplier that gives a linear increase in pore
%           volume with pressure, multiplied with cR::
%
%               pv = pv_ref * (1 + (p-pRef)*cR)
%
%           where pv_ref is the pore volume specified by the grid and rock
%           models
%
%   smin - Minimal (connate) saturations for each phase. If > 0, these will 
%          scale rel-perms such that krX=0 for s <= smin(X) and krX = 1 for
%          s >= 1 - sum(smin) + smin(X)
%
%   pRef - Reference pressure used in conjunction with rock and fluid
%          compressibility factors c and cr, i.e::
%
%               b(p) = b_ref * exp((p-p_ref)*c)
%               pv = pv_ref * (1 + (p-pRef)*cR)
%
%          Default value is `pRef = 0`
%
% RETURNS:
%   fluid - struct containing the following functions (where X = 'W' [water],
%           'O' [oil] and 'G' [gas]):
%             * krW, krO, krG - Relative permeability functions
%             * rhoXS         - density of X at surface conditions
%             * bX(p)         - inverse formation volume factor
%             * muX(p)        - viscosity functions (constant)
%             * krX(s)        - rel.perm for X
%             * krOW, krOG    - (If 3ph - oil-water and oil-gas rel.perm)
%
% EXAMPLE: 
%   % Create an incompressible three-phase fluid with properties:
%   % water: rhoS = 1000kg/m^3, mu = 1.0cP
%   % oil:   rhoS =  700kg/m^3, mu = 5.0cP
%   % gas :  rhoS =  200kg/m^3, mu = 0.1cP
%   % And linear relative permeabilities
%   fluid = initSimpleADIFluid('phases', 'WOG', ...
%                              'mu',      [1, 5, 0.1]*centi*poise,...
%                              'rho'     [1000, 700, 200]*kilogram/meter^3)
% 
%   % Create a two-phase, oil-water fluid where the oil phase is weakly
%   % compressible and quadratic relative permeabilities for both phases:
%   % water: rhoS = 1000kg/m^3, mu = 1.0cP, incompressible
%   % oil:   rhoS =  650kg/m^3, mu = 8.0cP, b(p)=exp((p-100*barsa)*1e-5/barsa)
%
%   fluid = initSimpleADIFluid('phases', 'OW', ...
%                              'mu',      [8, 1]*centi*poise, ...
%                              'rho',     [650, 1000]*kilogram/meter^3 , ...
%                              'n',       [2, 2], ...
%                              'pRef',    100*barsa, ...
%                              'c',       [1e-5, 0]/barsa);

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

opt = struct('mu',       [], ...
             'rho',      [1, 1, 1], ...
             'n',        [1, 1, 1], ...
             'b',        [1, 1, 1], ...
             'c',        [], ...
             'pRef',     0, ...
             'cR',       [], ...
             'blackoil', true, ...
             'smin',     [0 0 0], ...
             'phases',   'WOG');
opt = merge_options(opt, varargin{:});
if isempty(opt.mu) && opt.blackoil
    opt.mu = [1, 1, 1];
    warning('Defaulted viscosity to unit value (1000 centipoise) for all phases.')
end

names = upper(opt.phases);
assert(all(ismember(double(names), double('WOG'))), ...
    'Input phases must be some combination of W for water, O for oil and G for gas!');
nPh = numel(names);
assert(nPh == numel(unique(double(names))), 'Duplicate phases detected.');
assert(sum(opt.smin) < .99, 'Sum of connate/critical saturations must be smaller than 1');
for i = 1:nPh
    n = names(i);
    b = opt.b(i);
    if isempty(opt.c)
        % Constant value (incompressible phase)
        bf = @(p, varargin) b*constantReciprocalFVF(p, varargin{:});
    else
        % Compressibility on the form
        % b = b_ref exp((p-p_ref)*c)
        c = opt.c(i);
        if c < 0
            warning('Negative compressibility detected.')
        end
        bf = @(p, varargin) b*exp((p-opt.pRef)*c);
    end
    if all(opt.smin == 0)
        kr = @(s, varargin) s.^opt.n(i);
        [sl, su] = deal(0,1);
    else
        sl = opt.smin(i);
        su = 1 - sum(opt.smin) + sl;
        kr = @(s, varargin) max(0, min(1, (s-sl)/(su-sl) )).^opt.n(i);
    end
    
    fluid.(['rho', n, 'S']) = opt.rho(i);
    if opt.blackoil
        fluid.(['b', n]) = bf;
        fluid.(['mu', n]) = @(p, varargin) constantViscosity(opt.mu(i), p, varargin{:});
    end
    fluid.(['kr', n]) = kr;
    fluid.krPts.(lower(n)) = [sl, sl, su, 1];
    if strcmpi(n, 'O') && nPh > 2
        [fluid.krOW, fluid.krOG] = deal(kr);
        [fluid.krPts.ow, fluid.krPts.og] = deal([sl, sl, su, 1]);
    end
end

if ~isempty(opt.cR)
    % Rock compressibility
    cR = opt.cR;
    assert(numel(cR) == 1, 'Rock compressibility must be given as a single number');
    assert(cR >= 0, 'Rock compressibility must be a positive number');
    fluid.cR = cR;
    fluid.pvMultR = @(p)(1 + cR.*(p-opt.pRef));
end


end

function B = constantReciprocalFVF(p, varargin)
B = p.*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
mu = p.*0 + mu;
end
