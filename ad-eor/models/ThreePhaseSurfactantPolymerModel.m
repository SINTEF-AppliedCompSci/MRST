classdef ThreePhaseSurfactantPolymerModel < ThreePhaseBlackOilModel
    % Three-phase black-oil model with support for surfactant and polymer injection
    %
    % SYNOPSIS:
    %   model = ThreePhaseSurfactantPolymerModel(G, rock, fluid, varargin)
    %
    % DESCRIPTION:
    %   Fully implicit three phase blackoil model with polymer.
    %
    % PARAMETERS:
    %   G        - Grid
    %   rock     - Rock structure
    %   fluid    - Fluid structure
    %   varargin - optional parameters
    %
    % RETURNS:
    %   class instance
    %
    % EXAMPLE:
    %
    % SEE ALSO:  equationsThreePhaseBlackOilPolymer, OilWaterPolymerModel
    %

    properties
        % Surfactant and Polymer present
        polymer
        surfactant
        % Using PLYSHEAR shear model based on water velocity
        usingShear
	    % Using PLYSHLOG shear model based on water velocity
        usingShearLog
	    % Using PLYSHLOG shear model based on water shear rate
        usingShearLogshrate
    end

    methods
        function model = ThreePhaseSurfactantPolymerModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid, varargin{:});

            if isempty(model.inputdata)
                % We guess what's present
                model.polymer = isfield(model.fluid, 'cpmax');
                model.surfactant = isfield(model.fluid, 'ift');
            else
                % We have a deck that explicitly enables features
                runspec = model.inputdata.RUNSPEC;
                check = @(name) isfield(runspec, upper(name)) && runspec.(upper(name));
                model.polymer    = check('POLYMER');
                model.surfactant = check('SURFACT');
            end

            hasSHRATE    = isfield(fluid, 'shrate');
            hasPLYSHLOG  = isfield(fluid, 'plyshlog');
            hasPLYSHMULT = isfield(fluid, 'plyshearMult');

            if (hasSHRATE && ~hasPLYSHLOG)
                error('SHRATE is specified but PLYSHLOG is not')
            end
            if (hasPLYSHMULT && hasPLYSHLOG)
                error('PLYSHLOG and PLYSHEAR are existing together');
            end

            model.usingShear    = hasPLYSHMULT;
            model.usingShearLog = hasPLYSHLOG && ~hasSHRATE;
            model.usingShearLogshrate =  hasPLYSHLOG && hasSHRATE;
            model = merge_options(model, varargin{:});
        end

        % --------------------------------------------------------------------%
        function model = setupOperators(model, G, rock, varargin)
            model.operators = setupOperatorsTPFA(G, rock, varargin{:});
            if model.surfactant
                % Operators used to compute capillary number
                model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
                model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
            end
        end

        % --------------------------------------------------------------------%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            [problem, state] = equationsThreePhaseSurfactantPolymer(state0, state, ...
                model, dt, drivingForces, varargin{:});
        end
        % --------------------------------------------------------------------%
        function model = validateModel(model, varargin)
            model = validateModel@ThreePhaseBlackOilModel(model, varargin{:});
            if model.surfactant && ~isfield(model.operators, 'veloc')
                % Operators used to compute capillary number
                G = model.G;
                model.operators.veloc = computeVelocTPFA(G, model.operators.internalConn);
                model.operators.sqVeloc = computeSqVelocTPFA(G, model.operators.internalConn);
            end
        end
        
        % --------------------------------------------------------------------%
        function state = validateState(model, state)
            state = validateState@ThreePhaseBlackOilModel(model, state);
            nc = model.G.cells.num;
            if model.polymer
                model.checkProperty(state, 'Polymer', [nc, 1], [1, 2]);
                fn = model.getVariableField('Polymermax');
                if ~isfield(state, fn)
                    state.(fn) = model.getProp(state, 'Polymer');
                end
                model.checkProperty(state, 'Polymermax', [nc, 1], [1, 2]);
            end
            if model.surfactant
                model.checkProperty(state, 'Surfactant', [nc, 1], [1, 2]);
                fn = model.getVariableField('SurfactantMax');
                if ~isfield(state, fn)
                    state.(fn) = model.getProp(state, 'Surfactant');
                end
                model.checkProperty(state, 'SurfactantMax', [nc, 1], [1, 2]);
            end
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            [state, report] = updateState@ThreePhaseBlackOilModel(model, ...
               state, problem,  dx, drivingForces);

            % cp denotes concentration of polymer, cs is concentration of surfactant.
            if model.polymer
                cp = model.getProp(state, 'polymer');
                cp = min(cp, model.fluid.cpmax);
                state = model.setProp(state, 'polymer', max(cp, 0));
            end
            if model.surfactant
                cs = model.getProp(state, 'surfactant');
                state = model.setProp(state, 'surfactant', max(cs, 0) );
            end
        end

        % --------------------------------------------------------------------%
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            [state, report] = updateAfterConvergence@ThreePhaseBlackOilModel(model, state0, state, dt, drivingForces);

            if model.polymer
                cp    = model.getProp(state, 'polymer');
                cpmax = model.getProp(state, 'polymermax');
                state = model.setProp(state, 'polymermax', max(cpmax, cp));
            end
            if model.surfactant
                cs    = model.getProp(state, 'surfactant');
                csmax = model.getProp(state, 'surfactantmax');
                state = model.setProp(state, 'surfactantmax', max(csmax, cs));
            end
        end

        % --------------------------------------------------------------------%
        function model = setupStateFunctionGroupings(model, varargin)
            model = setupStateFunctionGroupings@ThreePhaseBlackOilModel(model, varargin{:});

            fp = model.FlowPropertyFunctions;
            pp = model.PVTPropertyFunctions;
            fd = model.FlowDiscretization;
            pvtreg  = pp.getRegionPVT(model);
            satreg  = fp.getRegionSaturation(model);
            surfreg = fp.getRegionSurfactant(model);

            % We set up EOR viscosities and relative permeabilities. They are computed from
            % the black-oil value by using a multiplier approach, where we have
            % one multiplier for each EOR effect.
            pp = pp.setStateFunction('Viscosity', EORViscosity(model, pvtreg));
            pp = pp.setStateFunction('BaseViscosity', BlackOilViscosity(model));
            fp = fp.setStateFunction('RelativePermeability', EORRelativePermeability(model));
            fp = fp.setStateFunction('BaseRelativePermeability', BaseRelativePermeability(model, satreg));

            % The statefunction ViscosityMultipliers and RelPermMultipliers are containers
            % for the Viscosity and Relative Permeability multpliers.  Each
            % multiplier is set up as a property (for example
            % 'PolymerEffViscMult' below) and added to the container.
            viscmult = PhaseMultipliers(model);
            viscmult.label = 'M_\mu';

            % TODO: we need to find an easier way to handle the viscosity multiplication between different component
            % hopefully, we only handle once.
            pviscmult = PhaseMultipliers(model);
            pviscmult.label = 'M_{\mu_p}';

            relpermult = PhaseMultipliers(model);
            relpermult.label = 'M_{kr}';
            relpermult.operator = @rdivide; % The relperm multipliers are divided

            if model.polymer

                fp = fp.setStateFunction('PolymerAdsorption', PolymerAdsorption(model, satreg));
                fd = fd.setStateFunction('PolymerPhaseFlux' , PolymerPhaseFlux(model));
                fd = fd.setStateFunction('FaceConcentration', FaceConcentration(model));
                fd = fd.setStateFunction('ComponentPhaseFlux', ComponentPhaseFluxWithPolymer(model));

                % We set up the water effective viscosity multiplier based on polymer concentration
                peffmult = 'PolymerEffViscMult';
                pp = pp.setStateFunction(peffmult, PolymerEffViscMult(model, pvtreg));
                viscmult = viscmult.addMultiplier(model, peffmult, 'W');

                % We set up the polymer effective viscosity multiplier, which is used for polymer transport
                polyviscmult = 'PolymerViscMult';
                pp = pp.setStateFunction(polyviscmult, PolymerViscMult(model, pvtreg));
                % TODO: really not sure whether this is the correct way to
                % do pviscmult
                pviscmult = pviscmult.addMultiplier(model, polyviscmult, 'W');

                % We set up the polymer permeability reduction effect. If permeability reduction
                % is present, it means that we must divide the water relative
                % permeability with the value of the "multiplier".
                rppmult = 'PolymerPermReduction';
                fp = fp.setStateFunction(rppmult, PolymerPermReduction(model));
                relpermult = relpermult.addMultiplier(model, rppmult, 'W');
                % if necessary, the following can be moved out similar with
                % others, currently they are overwritten due to polymer
                if ~isempty(model.FacilityModel) && isprop(model.FacilityModel, 'FacilityFlowDiscretization')
                    ffd = model.FacilityModel.FacilityFlowDiscretization;
                    ffd = ffd.setStateFunction('Mobility', PerforationMobilityEOR(model));
                    if model.polymer
                        ffd = ffd.setStateFunction('ComponentPhaseDensity', PerforationComponentPhaseDensityEOR(model));
                    end
                    model.FacilityModel.FacilityFlowDiscretization = ffd;
                end
            end

            if model.surfactant
                fp = fp.setStateFunction('CapillaryNumber', CapillaryNumber(model));
                fp = fp.setStateFunction('SurfactantAdsorption', SurfactantAdsorption(model, satreg));
                % The EOR relative permeability is set up as the
                % SurfactantRelativePermeability combined with multipliers.
                fp = fp.setStateFunction('BaseRelativePermeability', SurfactantRelativePermeability(model, satreg, surfreg));
                fp.CapillaryPressure = SurfactantCapillaryPressure(model, satreg);

                % We set up the surfactant viscosity multiplier
                smult = 'SurfactantViscMultiplier';
                pp = pp.setStateFunction(smult, SurfactantViscMultiplier(model, pvtreg));
                viscmult = viscmult.addMultiplier(model, smult, 'W');
                pviscmult = pviscmult.addMultiplier(model, smult, 'W');
            end

            pp = pp.setStateFunction('ViscosityMultipliers', viscmult);
            % TODO: BAD NAME!
            pp = pp.setStateFunction('PolyViscMult', pviscmult);
            fp = fp.setStateFunction('RelativePermeabilityMultipliers', relpermult);

            model.FlowPropertyFunctions = fp;
            model.PVTPropertyFunctions  = pp;
            model.FlowDiscretization    = fd;            
        end

        % --------------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            % Get the index/name mapping for the model (such as where
            % pressure or water saturation is located in state)
            index = ':';
            switch(lower(name))
                case {'polymer', 'cp'}
                    fn = 'cp';
                case {'polymermax', 'cpmax'}
                    fn = 'cpmax';
                case {'surfactant', 'cs'}
                    fn = 'cs';
                case {'surfactantmax', 'csmax'}
                    fn = 'csmax';
                case 'qwpoly'
                    fn = 'qWPoly';
                case 'qwsft'
                    fn = 'qWSft';
                otherwise
                    [fn, index] = getVariableField@ThreePhaseBlackOilModel(...
                                    model, name, varargin{:});
            end
        end

        % --------------------------------------------------------------------%
        function names = getComponentNames(model)
            names = getComponentNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'polymer';
            end
            if model.surfactant
                names{end+1} = 'surfactant';
            end
        end

        % --------------------------------------------------------------------%
        function state = storeSurfData(model, state, s, cs, Nc, sigma)
            state.SWAT    = double(s);
            state.SURFACT = double(cs);
            state.SURFCNM = log(double(Nc))/log(10);
            state.SURFST  = double(sigma);
            % state.SURFADS = double(ads);
        end

        % --------------------------------------------------------------------%
        function scaling = getScalingFactorsCPR(model, problem, names, solver)
            nNames = numel(names);

            scaling = cell(nNames, 1);
            handled = false(nNames, 1);

            for iter = 1:nNames
                name = lower(names{iter});
                switch name
                    case 'polymer'
                        s = 0;
                    otherwise
                        continue
                end
                sub = strcmpi(problem.equationNames, name);

                scaling{iter} = s;
                handled(sub) = true;
            end
            if ~all(handled)
                % Get rest of scaling factors
                other = getScalingFactorsCPR@ThreePhaseBlackOilModel(model, problem, names(~handled), solver);
                [scaling{~handled}] = other{:};
            end
        end

        % --------------------------------------------------------------------%
        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % PARAMETERS:
        %
        %   model  - (Base class, automatic)
        %
        %   cname  - Name of the component. Must be a property known to the
        %            model itself through `getProp` and `getVariableField`.
        %
        %   eq     - Equation where the source terms are to be added. Should
        %            be one value per cell in the simulation grid (model.G)
        %            so that the src.sourceCells is meaningful.
        %
        %   component - Cell-wise values of the component in question. Used
        %               for outflow source terms only.
        %
        %   src    - Source struct containing fields for fluxes etc. Should
        %            be constructed from force and the current reservoir
        %            state by `computeSourcesAndBoundaryConditionsAD`.
        %
        %   force  - Force struct used to produce src. Should contain the
        %            field defining the component in question, so that the
        %            inflow of the component through the boundary condition
        %            or source terms can accurately by estimated.
            if isempty(force)
                return
            end
            c = model.getProp(force, cname);
            cells = src.sourceCells;
            switch lower(cname)
              case {'polymer'}
                % Water based EOR, multiply by water flux divided by
                % density and add into corresponding equation
                qW = src.phaseMass{1}./model.fluid.rhoWS;
                isInj = qW > 0;
                qC = (isInj.*c + ~isInj.*component(cells)).*qW;
              otherwise
                error(['Unknown component ''', cname, '''. BC not implemented.']);
            end
            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end

        function [names, types] = getExtraWellEquationNames(model)
            [names, types] = getExtraWellEquationNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'polymerWells';
                types{end+1} = 'perf';
            end
            if model.surfactant
                names{end+1} = 'surfactantWells';
                types{end+1} = 'perf';
            end
        end

        function names = getExtraWellPrimaryVariableNames(model)
            names = getExtraWellPrimaryVariableNames@ThreePhaseBlackOilModel(model);
            if model.polymer
                names{end+1} = 'qWPoly';
            end
            if model.surfactant
                names{end+1} = 'qWSft';
            end
        end

        function [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration)
            [compEqs, compSrc, eqNames, wellSol] = getExtraWellContributions@ThreePhaseBlackOilModel(model, well, wellSol0, wellSol, q_s, bh, packed, qMass, qVol, dt, iteration);
            if model.polymer || model.surfactant
                assert(model.water, 'Surfactant-Polymer injection requires a water phase.');
                f = model.fluid;

                % Water is always first
                wix = 1;
                cqWs = qMass{wix}./f.rhoWS; % connection volume flux at surface condition
            end

            if model.polymer
                if well.isInjector()
                    concWellp = model.getProp(well.W, 'polymer');
                    cqP = concWellp.*cqWs;
                else
                    pixp = strcmpi(model.getComponentNames(), 'polymer');
                    concWellp = packed.components{pixp};

                    a = f.muWMult(f.cpmax).^(1-f.mixPar);
                    cpbarw = concWellp/f.cpmax;

                    % the term (a + (1 - a).*cpbarw) account for the
                    % todd-longstaff mixing factor, which model the fact that for
                    % not-fully mixed polymer solution the polymer does not
                    % travel at the same velocity as water. See the governing
                    % equation for polymer (e.g. equationsOilWaterPolymer.m)
                    cqP = concWellp.*cqWs./(a + (1-a).*cpbarw);
                end
                qwpoly = packed.extravars{strcmpi(packed.extravars_names, 'qwpoly')};

                compEqs{end+1} = qwpoly - sum(concWellp.*cqWs);
                compSrc{end+1} = cqP;
                eqNames{end+1} = 'polymerWells';
            end

            if model.surfactant
                if well.isInjector()
                    concWells = model.getProp(well.W, 'surfactant');
                else
                    pixs = strcmpi(model.getComponentNames(), 'surfactant');
                    concWells = packed.components{pixs};
                end
                cqS = concWells.*cqWs;
                qwsft = packed.extravars{strcmpi(packed.extravars_names, 'qwsft')};

                compEqs{end+1} = qwsft - sum(concWells.*cqWs);
                compSrc{end+1} = cqS;
                eqNames{end+1} = 'surfactantWells';
            end
        end
    end
end

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
