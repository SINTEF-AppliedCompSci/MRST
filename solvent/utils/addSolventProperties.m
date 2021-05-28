function fluid = addSolventProperties(fluid, varargin)
% Add solvent pseudo-component and properties to fluid.

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

opt = struct('muS'             , 1    , ... % Solvent viscosity
             'rhoSS'           , 1    , ... % Solvent surface density
             'bS'              , 1    , ... % Formation volume factor
             'c'               , []   , ... % Compressibility
             'pRef'            , 0    , ... % Reference pressure
             'mixPar'          , 1    , ... % Mixing parameter for viscosity
             'mixParRho'       , []   , ... % Mixing parameter for density
             'sWcon'           , 0    , ... % Residual water saturation
             'sOr_m'           , 0    , ... % Miscible residual oil saturation
             'sGc_m'           , 0    , ... % Miscible residual gas saturation
             'Ms'              , []   , ... % Saturation-dependent miscibility
             'Mp'              , []   , ... % Pressure-dependent miscibility
             'MkrO'            , []   , ... % Miscible oil relperm multiplier function
             'MkrG'            , []   , ... % Miscible gas relperm multiplier function
             'krFS'            , []   , ... % Immiscible gas relperm multiplier function
             'krFG'            , []   , ... % Miscible solvent relperm multiplier function
             'smin'            , 1e-13, ... % Cut-off used to avoid division by zero
             'overwrite'       , false);    % Overwrite any solvent propertis already defined

opt = merge_options(opt, varargin{:});

%% Set solvent four-phase solvent model specifics

    overwrite = opt.overwrite;
    
    if isempty(opt.mixParRho)
        opt.mixParRho = opt.mixPar;
    end

    % Mixing parameters for viscosity and density
    
    fluid.sOr_i = getResSat(fluid.krOW);
    fluid.sGc_i = getResSat(fluid.krG);
    if overwrite || ~isfield(fluid, 'sWcon')
        fluid.sWcon = opt.sWcon;
    end
    
    if isempty(opt.Mp)
        opt.Mp = 1;
    end
    
    names = fieldnames(opt);
    for fNo = 1:numel(names)
        name = names{fNo};
        if any(strcmpi(name, {'mixPar', 'mixParRho', 'smin'}))
            type = {'constant'};
        elseif any(strcmpi(name, {'sOr_m', 'sGc_m', 'Mp'}))
            slope = 0; x0 = opt.(name);
            type = {'linear', slope, x0};
        elseif any(strcmpi(name, {'Ms', 'MkrG', 'krFS', 'krFG'}))
            slope = 1; x0 = 0;
            type = {'linear', slope, x0};
        elseif strcmp(name, 'MkrO')
            slope = -1; x0 = 1;
            type = {'linear', slope, x0};
        else
            continue
        end
        fluid = assignProp(fluid, names{fNo}, opt.(names{fNo}), type, overwrite);
    end
        
    fluid.satFrac = @(sX, sY) satFrac(sX, sY, fluid.smin);

%% Set standard properties for solvent "pahse" that are not already defined

if overwrite || ~isfield(fluid, 'bS')
    b = opt.bS;
    if isempty(opt.c)
        % Constant value (incompressible phase)
        bS = @(p, varargin) b*constantReciprocalFVF(p, varargin{:});
    else
        % Compressibility on the form
        % b = b_ref exp((p-p_ref)*c)
        c = opt.c;
        if c < 0
            warning('Negative compressibility detected.')
        end
        bS = @(p, varargin) b*exp((p-opt.pRef)*c);
    end
    fluid.bS    = bS;
end

if overwrite || ~isfield(fluid, 'rhoSS')
    fluid.rhoSS = opt.rhoSS;
end

if overwrite || ~isfield(fluid, 'muS')
    fluid.muS   = @(p, varargin) constantViscosity(opt.muS, p, varargin{:});
end

end

function fluid = assignProp(fluid, name, value, type, overwrite)
    if overwrite || ~isfield(fluid, name)
        fluid.(name) = value;
        switch type{1}
            case 'linear'
                if ~isa(fluid.(name), 'function_handle')
                    fluid.(name) = setLinearFunction(type{2}, type{3});
                end
        end
    end
end

function F = satFrac(sX, sY, smin)
    F     = sX;
    ii    = sX + sY > smin;
    if isnumeric(F) && isa(sY, 'ADI')
        if isa(sY, 'GenericAD')
            F = double2GenericAD(F, sY);
        else
            F = double2ADI(F, sY);
        end
    end
    F(ii) = sX(ii)./sY(ii);
    F = min(max(F, 0),1);
    
end

function B = constantReciprocalFVF(p, varargin)
    B = p*0 + 1;
end

function mu = constantViscosity(mu, p, varargin)
    mu = p*0 + mu;
end

function v = setLinearFunction(slope, x0)
    v = @(x) min(max(slope*x + x0,0),1);
end

function sr = getResSat(kr)

    s = linspace(0,1,101)';
    kr = kr(s);
    is_mob = kr > 0;
    sr = s([is_mob(2:end); true] & ~is_mob);
    
    if isempty(sr) 
        sr = 0;
    end
    
end