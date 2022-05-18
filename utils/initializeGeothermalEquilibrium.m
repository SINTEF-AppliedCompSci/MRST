function [p,T,sol] = initializeGeothermalEquilibrium(model, varargin)
%Compute geothermal (hydriostatic and thermal) equilibrium using the
%explicit Runge-Kutta (4,5) solver provided with ode45.
%
% SYNOPSIS:
%       [p,T,sol] = initializeGeothermalEquilibrium(model)
%       [p,T,sol] = initializeGeothermalEquilibrium(model, 'pn1', pv1, ...)
%
% PARAMETERS:
%   model - MRST model defining the fluid physics needed for computing the
%           equilibrium
%
% OPTIONAL ARGUMENTS:
%   'datumDepth'    - Depth at which the datum pressure and temperature is
%                  given. Default is 0 (surface level).
%
%   'datumPressure/Temperature' - Pressure/temperature at the datum depth.
%   Detault is 1*atm/(273.15 + 20)*Kelvin
%
%   'thermalGradient' - Thermal gradient as a function of depth. If left
%   empty, the function uses model.thermalGradient if this exists, and
%   30*Kelvin/(kilo*meter) if not.
%
% RETURNS:
%   p,T - Equilibrium pressure and temperature for each grid cell
%
%   sol - The solution for p/T as returned from ode45. The solution can be
%         evaluated in z coordinates as [p,T] = sol(z)
%
% SEE ALSO:
%   `GeothermalModel`

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

    % Optional input arguments
    %---------------------------------------------------------------------%
    opt = struct('datumDepth'      , 0*meter           , ...
                 'datumPressure'   , 1*atm             , ...
                 'datumTemperature', 273.15 + 20*Kelvin, ...
                 'thermalGradient' , []                );
    opt = processInput(model, opt, varargin{:});
    %---------------------------------------------------------------------%
    
    % Set up and solve ODE defining geothermal equilibrilum
    %---------------------------------------------------------------------%
    [dpdz, dTdz, p0, T0, z, z0] = getPhysicalQuantities(model, opt);
    [p,T,sol] = solveODE(dpdz, dTdz, p0, T0, z, z0);
    %---------------------------------------------------------------------%
end

%-------------------------------------------------------------------------%
function opt = processInput(model, opt, varargin)
% Process input arguments

    opt = merge_options(opt, varargin{:});
    if isempty(opt.thermalGradient)
        if isprop(model, 'thermalGradient')
            opt.thermalGradient = model.thermalGradient;
        else
            opt.thermalGradient = 30*Kelvin/(kilo*meter);
        end
    end
    if isscalar(opt.thermalGradient)
        opt.thermalGradient = @(p,T) opt.thermalGradient;
    end
    assert(isa(opt.thermalGradient, 'function_handle'), ...
                    'Thermal gradient must be a scalar or function handle');

end

%-------------------------------------------------------------------------%
function [dpdz, dTdz, p0, T0, z, z0] = getPhysicalQuantities(model, opt)
% Get physical quantities defining the ODE

    % Pressure gradient due to density
    g    = norm(model.gravity);
    rho  = @(p,T) model.fluid.rhoW(p, T);
    dpdz = @(p,T) g*rho(p,T);
    % Tempreature gradient if constant
    dTdz = opt.thermalGradient();
    % Datum pressure and temperature
    p0 = opt.datumPressure;
    T0 = opt.datumTemperature;
    % Verical coordinates in the grid
    z  = model.G.cells.centroids(:,3);
    z0 = opt.datumDepth;

end

%-------------------------------------------------------------------------%
function [p,T,sol] = solveODE(dpdz, dTdz, p0, T0, z, z0)
% Solve ODE to find equilibrium pressure and temperature

    % Set initial condition
    y0 = [p0; T0];
    % Set derivative
    dydz = @(z,y) [dpdz(y(1), y(2)); dTdz + y(2)*0];
    % Set ode45 options
    odeopts = odeset('AbsTol', 1.0e-10, 'RelTol', 5.0e-8);
    % Solve ODE
    [solUp, solDown] = deal([]);
    if z0 > min(z)
        zlimUp = [z0, min(z)];
        solUp  = ode45(dydz, zlimUp, y0, odeopts);
    end
    if z0 < max(z)
        zlimDown = [z0, max(z)];
        solDown  = ode45(dydz, zlimDown, y0, odeopts);
    end
    assert(~(isempty(solUp) && isempty(solDown)), 'Something went wrong');
    % Evaluate in cells
    [p,T] = evaluateFn(solUp, solDown, z0, z);
    sol   = @(z) evaluateFn(solUp, solDown, z);
    
end

%-------------------------------------------------------------------------%
function [p, T] = evaluateFn(solUp, solDown, z0, z)
% Evaluate solution in given z coordinates

    [p,T] = deal(zeros(size(z)));
    
    top = z <= z0; % Include endpoint in case z0 == max(z)
    if any(top) && ~isempty(solUp)
        y = deval(solUp, z(top)); p(top) = y(1,:); T(top) = y(2,:);
    end
    
    bot = z >= z0; % Include endpoint in case z0 == min(z)
    if any(bot) && ~isempty(solDown)
        y = deval(solDown, z(bot)); p(bot) = y(1,:); T(bot) = y(2,:)';
    end
    
end