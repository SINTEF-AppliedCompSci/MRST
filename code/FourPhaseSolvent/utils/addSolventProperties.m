function fluid = addSolventProperties(fluid, varargin)
% Add solvent "phase" and properties to fluid.

%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

opt = struct('mu'              , 1    , ... % Solvent viscosity
             'rho'             , 1    , ... % Solvent surface density
             'n'               , 1    , ... % Solvent Relperm coefficient
             'b'               , 1    , ... % Formation volume factor
             'c'               , []   , ... % Compressibility
             'pRef'            , 0    , ... % Reference pressure
             'mixPar'          , 1    , ... % Mixing parameter
             'sWr'             , 0    , ... % Residual water saturation
             'sOr_m'           , 0    , ... % Miscible residual oil saturation
             'sOr_i'           , 0    , ... % Immiscible residual oil saturation
             'sGc_m'           , 0    , ... % Miscible residual gas saturation
             'sGc_i'           , 0    , ... % Immiscible residual gas saturation
             'Ms'              , []   , ... % Saturation-dependent miscibility
             'Mp'              , []   , ... % Pressure-dependent miscibility
             'MkrO'            , []   , ... % Miscible oil relperm multiplier function
             'MkrG'            , []   , ... % Miscible gas relperm multiplier function
             'krFS'            , []   , ... % Immiscible gas relperm multiplier function
             'krFG'            , []   , ... % Miscible solvent relperm multiplier function
             'smin'            , 1e-10);    % Cut-off used to avoid division by zero

opt = merge_options(opt, varargin{:});

%% Set solvent four-phase solvent model specifics

fluid.mixPar = opt.mixPar;
fluid.smin  = opt.smin;
fluid.sWr    = opt.sWr;

% Miscible and immiscible residual oil and gas saturations. Miscible
% residual saturation may be an increasing function of water saturation.
fluid.sOr_i = opt.sOr_i;
fluid.sOr_m = opt.sOr_m;
if ~isa(fluid.sOr_m, 'functionHandle')
    fluid.sOr_m = setConstantFunction(fluid.sOr_m);
end

fluid.sGc_i  = opt.sGc_i;
fluid.sGc_m = opt.sGc_m;
if ~isa(fluid.sGc_m, 'functionHandle')
    fluid.sGc_m = setConstantFunction(fluid.sGc_m);
end

% Saturation-dependent miscibiliy function. 0 <= Ms(sS/(sS + sG) <= 1,
% increasing.
fluid.Ms = opt.Ms;
if isempty(fluid.Ms)
    fluid.Ms = setLinearFunction();
end

% Saturation-dependent miscibiliy function. 0 <= Mp(p) <= 1, increasing.
fluid.Mp = opt.Mp;
if isempty(fluid.Mp)
    fluid.Mp = setConstantFunction(1);
end

% Miscible relperm multipliers.
% krO_i = MkrO((sO-sOr)/(sO + sG + sS -(sOr + sGc)))*krN(sO + sG + sS) etc
fluid.MkrO = opt.MkrO;
if isempty(fluid.MkrO)
    fluid.MkrO = setLinearFunction();
end

fluid.MkrG = opt.MkrG;
if isempty(fluid.MkrG)
    fluid.MkrG = setLinearFunction();
end

% Immiscible relperm multipliers.
% krG_m = krFG((sG-sGc)/(sG + sS - sGc))*krG(sG + sS) etc
fluid.krFS = opt.krFS;
if isempty(fluid.krFS)
    fluid.krFS = setLinearFunction();
end

fluid.krFG = opt.krFG;
if isempty(fluid.krFG)
    fluid.krFG = setLinearFunction();
end

% Model equaitons involves expressions of the form sA/(sA + sB), which may
% result in singular expressions. Function satFrac gracefully handels this
% for us.
fluid.satFrac = @(sX, sY) satFrac(sX, sY, fluid.smin);

%% Set standard properties (b, kr, rho, mu) for the solvent "pahse"

b = opt.b;
if isempty(opt.c)
    % Constant value (incompressible phase)
    bf = @(p, varargin) b*constantReciprocalFVF(p, varargin{:});
else
    % Compressibility on the form
    % b = b_ref exp((p-p_ref)*c)
    c = opt.c;
    if c < 0
        warning('Negative compressibility detected.')
    end
    bf = @(p, varargin) b*exp((p-opt.pRef)*c);
end
    
fluid.rhoSS = opt.rho;
fluid.bS    = bf;
fluid.muS   = @(p, varargin) constantViscosity(opt.mu, p, varargin{:});


end

function F = satFrac(sX, sY, smin)
    F     = sX;
    is    = sX + sY > smin;
    F(is) = sX(is)./sY(is);
end

function B = constantReciprocalFVF(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end

function v = setConstantFunction(v)
    v = @(x) x.*0 + v;
end

function v = setLinearFunction()
    v = @(x) x;
end
    