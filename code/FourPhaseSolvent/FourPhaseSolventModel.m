classdef FourPhaseSolventModel < ReservoirModel
% Three phase with optional dissolved gas and vaporized oil
properties
   solvent
end

methods
    function model = FourPhaseSolventModel(G, rock, fluid, varargin)
        model = model@ReservoirModel(G, rock, fluid, varargin{:});

        % Blackoil -> use CNV style convergence 
        model.useCNVConvergence = true;

        % All phases are present
        model.water = true;
        model.oil = true;
        model.gas = true;
        model.solvent = true;
        model.saturationVarNames = {'sw', 'so', 'sg', 'ss'};

        model = merge_options(model, varargin{:});

    end
    
    function model = validateModel(model, varargin)
        if isempty(model.FacilityModel)
            model.FacilityModel = FacilityModelSolvent(model); %#ok
        end
        if nargin > 1
            W = varargin{1}.W;
            model.FacilityModel = model.FacilityModel.setupWells(W);
        end
        model = validateModel@PhysicalModel(model, varargin{:});
        return
    end
    
    % --------------------------------------------------------------------%
    function [fn, index] = getVariableField(model, name)
        switch(lower(name))
            case {'solvent', 'ss'}
                index = 4;
                fn = 's';
            otherwise
                % Basic phases are known to the base class
                [fn, index] = getVariableField@ReservoirModel(model, name);
        end
    end
    
    function vars = getSaturationVarNames(model)
        vars = {'sw', 'so', 'sg', 'ss'};
        ph = model.getActivePhases();
        vars = vars(ph);
    end
    
    % --------------------------------------------------------------------%
    function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        [problem, state] = equationsFourPhaseSolvent(state0, state, ...
                model, dt, drivingForces, varargin{:});

    end

    % --------------------------------------------------------------------%
    function state = validateState(model, state)
        % Check parent class
        state = validateState@ReservoirModel(model, state);
    end
    
    function phNames = getPhaseNames(model)
        % Get the active phases in canonical ordering
        tmp = 'WOGS';
        active = model.getActivePhases();
        phNames = tmp(active);
    end

    % --------------------------------------------------------------------%
    function [state, report] = updateState(model, state, problem, dx, drivingForces)
        saturations = lower(model.saturationVarNames);
        wi = strcmpi(saturations, 'sw');
        oi = strcmpi(saturations, 'so');
        gi = strcmpi(saturations, 'sg');
        si = strcmpi(saturations, 'ss');
        
        % Parent class handles almost everything for us
        [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);

        % Handle the directly assigned values (i.e. can be deduced directly from
        % the well controls. This is black oil specific.
        W = drivingForces.W;
        state.wellSol = assignWellValuesFromControlSolvent(model, state.wellSol, W, wi, oi, gi, si);
    end
    
    function [isActive, phInd] = getActivePhases(model)
        % Get active flag for the canonical phase ordering (water, oil
        % gas as on/off flags).
        isActive = [model.water, model.oil, model.gas, model.solvent];
        if nargout > 1
            phInd = find(isActive);
        end
    end
    
    function rhoS = getSurfaceDensities(model)
        active = model.getActivePhases();
        props = {'rhoWS', 'rhoOS', 'rhoGS', 'rhoSS'};
        rhoS = cellfun(@(x) model.fluid.(x), props(active));
    end
    
    function state = storeFluxes(model, state, vW, vO, vG, vS)
        % Utility function for storing the interface fluxes in the state
        isActive = model.getActivePhases();

        internal = model.operators.internalConn;
        state.flux = zeros(numel(internal), sum(isActive));
        phasefluxes = {double(vW), double(vO), double(vG), double(vS)};
        state = model.setPhaseData(state, phasefluxes, 'flux', internal);
    end
    
     % --------------------------------------------------------------------%
    function state = storeMobilities(model, state, mobW, mobO, mobG, mobS)
        % Utility function for storing the mobilities in the state
        isActive = model.getActivePhases();

        state.mob = zeros(model.G.cells.num, sum(isActive));
        mob = {double(mobW), double(mobO), double(mobG), double(mobS)};
        state = model.setPhaseData(state, mob, 'mob');
    end
    
    % --------------------------------------------------------------------%
    function state = storeUpstreamIndices(model, state, upcw, upco, upcg, upcs)
        % Store upstream indices, so that they can be reused for other
        % purposes.
        isActive = model.getActivePhases();

        nInterfaces = size(model.operators.N, 1);
        state.upstreamFlag = false(nInterfaces, sum(isActive));
        mob = {upcw, upco, upcg, upcs};
        state = model.setPhaseData(state, mob, 'upstreamFlag');
    end
    
    % --------------------------------------------------------------------%
    function state = storeDensity(model, state, rhoW, rhoO, rhoG, rhoS)
        % Store compressibility / surface factors for plotting and
        % output.
        isActive = model.getActivePhases();

        state.rho = zeros(model.G.cells.num, sum(isActive));
        rho = {double(rhoW), double(rhoO), double(rhoG), double(rhoS)};
        state = model.setPhaseData(state, rho, 'rho');
    end
    % --------------------------------------------------------------------%
    function state = storebfactors(model, state, bW, bO, bG, bS)
        % Store compressibility / surface factors for plotting and
        % output.
        isActive = model.getActivePhases();

        state.bfactor = zeros(model.G.cells.num, sum(isActive));
        b = {double(bW), double(bO), double(bG), double(bS)};
        state = model.setPhaseData(state, b, 'bfactor');
    end
    
end
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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