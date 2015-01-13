function fluid = initSimpleADIFluid(varargin)
% Make a structure representing a three-component fluid (water, oil, gas)
% and their properties (relative permeabilities, densities, viscosities).
% Relative permeabilities are modeled as monomial functions of saturation.
% Densities are assumed constant, so the returned formation volume factors
% are constant equal to 1.
%
% SYNOPSIS:
%   fluid = initSimpleADIFluid()
%   fluid = initSimpleADIFluid('pn1', 'pv1')
%
% PARAMETERS:
%   'pn'/pv - List of property names/property values.  Possibilities are:
%   - mu  : vector of viscosity values for water, oil and gas, [muW, muO, muG].  
%           Default is [1 1 1].
%   - rho : vector of density values for water, oil and gas, [rhoW, rhoO, rhoG].
%           Default is [1 1 1].
%   - n   : vector of the degrees of the monomials describing relative
%           permeability for water, oil and gas [nW, nO, nG].  Default is
%           [1 1 1] (linear relative permeabilities)
%
% RETURNS:
%   fluid - struct containing the following functions (where X = 'W' [water],
%           'O' [oil] and 'G' [gas]) 
%           * [krW, krO, krG] = relPerm(s_water, s_gas) - relative permeability functions         
%           * rhoX         - density of X
%           * rhoXS        - density of X at surface (equal to rhoX, since incompressible)
%           * bX(p), BX(p) - formation volume factors and their inverses
%                            (constants, always equal one)
%           * muX(p)       - viscosity functions (constant)
%           * krX(s)       - rel.perm for X
%           * krOX(s)      - 
%           * rsSat()      - saturation value for dissolved gas (always returns 0)

   opt = struct('mu', [1 1 1], 'rho', [1 1 1], 'n', [1 1 1]);
   opt = merge_options(opt, varargin{:});

   krW = @(sw, varargin) sw.^opt.n(1);
   krO = @(so, varargin) so.^opt.n(2);
   krG = @(sg, varargin) sg.^opt.n(3);
   relperms = {krW, krO, krG};

   fluid.relPerm = @(sw, sg, varargin) relPerm(krW, krO, krG, sw, sg, varargin{:});


   names = {'W', 'O', 'G'};
   for i = 1:numel(names)
       n = names{i};
       bf = @(p, varargin) constantUnitBfactor(p, varargin{:});

       fluid.(['rho', n]) = opt.rho(i);
       fluid.(['rho', n, 'S']) = opt.rho(i);
       fluid.(['b', n]) = bf;
       fluid.(['B', n]) = bf;
       fluid.(['mu', n]) = @(p, varargin) constantViscosity(opt.mu(i), p, varargin{:});
       fluid.(['kr', n]) = relperms{i};
       fluid.(['krO', n]) = relperms{i};
   end
   fluid.rsSat = @(varargin) varargin{1}*0;
end

function [krW, krO, krG] = relPerm(krW, krO, krG, sw, sg, varargin)
    krW = krW(sw, varargin{:});
    krO = krO(1 - sw - sg, varargin{:});
    krG = krG(sg, varargin{:});
end

function B = constantUnitBfactor(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end
