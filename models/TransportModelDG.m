classdef TransportModelDG < TransportModel
    
    properties
    end
    
    methods
        %-----------------------------------------------------------------%
        function model = TransportModelDG(parent, varargin)
           
            model = model@TransportModel(parent);
            
            [model, discArgs] = merge_options(model, varargin{:});
            % Construct discretization
            if isempty(model.disc)
                model.disc = DGDiscretization(model, discArgs{:});
            end
            model.parentModel.disc = model.disc;
        end
        
        %-----------------------------------------------------------------%
        function id = isDof(model, name)
            id = numel(name) > 3 && strcmp(name(end-2:end), 'dof');
        end
        
        %-----------------------------------------------------------------%
        function [fn, index] = getVariableField(model, name, varargin)
            [fn, index] = getVariableField@TransportModel(model, name, varargin{:});
            if ~isempty(fn)
                fn = [fn, 'dof'];
            end
        end
        
        %-----------------------------------------------------------------%
        function state = validateState(model, state)
            state    = assignDofFromState(model.disc, state);    
            state    = validateState@TransportModel(model, state);
            state = rmfield(state, 'sTdof');
            state.sT = sum(state.s,2);
            state    = assignDofFromState(model.disc, state);
        end
        
        function [state, names, origin] = getStateAD(model, state, init)
            if nargin < 3
                init = true;
            end
            parent = model.parentModel;
            % Get the AD state for this model
            [~, basenames, origin] = model.getPrimaryVariables(state);
            isParent = strcmp(origin, class(parent));
            basenames = basenames(isParent);
            origin = origin(isParent);
            basevars = cell(1, numel(basenames));
            for bNo = 1:numel(basenames)
                basevars{bNo} = model.getProp(state, basenames{bNo}, false);
            end
            % Find saturations
            isS = false(size(basevars));
            nph = parent.getNumberOfPhases();
            phase_variable_index = zeros(nph, 1);
            for i = 1:numel(basevars)
                [f, ix] = model.getVariableField(basenames{i});
                if strcmp(f, 'sdof')
                    isS(i) = true;
                    phase_variable_index(ix) = i;
                end
            end
            % Figure out saturation logic
            isP = strcmp(basenames, 'pressure');
            vars = basevars;
            names = basenames;
            useTotalSaturation = strcmpi(model.formulation, 'totalSaturation') ...
                                    && sum(isS) == nph - 1;
            if useTotalSaturation
                % Replace pressure with total saturation
                replacement = 'sT';
                sTdof = state.sTdof;
%                 sT = model.getProp(state, [replacement, 'dof'], false);
                % Replacing
                vars{isP} = sTdof;
                names{isP} = replacement;
                origin{isP} = class(model);
            else
                % Remove pressure and skip saturation closure
                vars = vars(~isP);
                names = names(~isP);
                origin = origin(~isP);
            end
            if init
                [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            end
            if useTotalSaturation
                basevars(~isP) = vars(~isP);
            else
                basevars(~isP) = vars;
            end
            state = model.initStateAD(state, basevars, basenames, origin);
            state.sdof = state.s;
            state.s = model.disc.getCellMean(state, value(state.sdof));
            if useTotalSaturation
                % Set total saturation as well
                sT = vars{isP};
                state = model.setProp(state, replacement, sT);
            end
        end
        
        %-----------------------------------------------------------------%
        function [eqs, names, types, state] = getModelEquations(model, state0, state, dt, drivingForces)
            tmodel = model.parentModel;
            [eqs, flux, names, types] = tmodel.FluxDiscretization.componentConservationEquations(tmodel, state, state0, dt);
            src = model.FacilityModel.getComponentSources(state);
            % Assemble equations and add in sources
            if strcmpi(model.formulation, 'missingPhase')
                % Skip the last phase! Only mass-conservative for
                % incompressible problems
                eqs = eqs(1:end-1);
                flux = flux(1:end-1);
                names = names(1:end-1);
                types = types(1:end-1);
            end
            for i = 1:numel(eqs)
                if ~isempty(src.cells)
                    eqs{i}(src.cells) = eqs{i}(src.cells) - src.value{i};
                end
                eqs{i} = model.operators.AccDiv(eqs{i}, flux{i});
                 if ~model.useCNVConvergence
                    pv     = model.operators.pv;
                    eqs{i} = eqs{i}.*(dt./pv);
                 end    
            end
        end
        
%         %-----------------------------------------------------------------%
%         function state = initStateAD(model, state, vars, names, origin)
%             removed = false(size(vars));
%             if model.disgas || model.vapoil
%                 % Black-oil specific variable switching
%                 if model.water
%                     isw = strcmpi(names, 'sw');
%                     sW = vars{isw};
%                     removed = removed | isw;
%                 else
%                     sW = 0;
%                 end
% 
%                 isx = strcmpi(names, 'x');
%                 x = vars{isx};
%                 sG = model.getProps(state, 'sg');
%                 st  = model.getCellStatusVO(state, 1-sW-sG, sW, sG);
%                 sG = st{2}.*(1-sW) + st{3}.*x;
%                 sO = st{1}.*(1-sW) + ~st{1}.*(1 - sW - sG);
%                 if model.water
%                     sat = {sW, sO, sG};
%                 else
%                     sat = {sO, sG};
%                 end
%                 removed(isx) = true;
%             else
%                 % Without variable switching
%                 phases = model.getPhaseNames();
%                 nph = numel(phases);
%                 sat = cell(1, nph);
%                 fill = ones(model.G.cells.num, 1);
%                 removed_sat = false(1, nph);
%                 for i = 1:numel(phases)
%                     sub = strcmpi(names, ['s', phases(i)]);
%                     if any(sub)
%                         fill = fill - vars{sub};
%                         removed = removed | sub;
%                         removed_sat(i) = true;
%                         sat{i} = vars{sub};
%                     end
%                 end
%                 if any(~removed_sat)
%                     sat{~removed_sat} = fill;
%                 end
%             end
%             state = model.setProp(state, 's', sat);
% 
%             if not(isempty(model.FacilityModel))
%                 % Select facility model variables and pass them off to attached
%                 % class.
%                 fm = class(model.FacilityModel);
%                 isF = strcmp(origin, fm);
%                 state = model.FacilityModel.initStateAD(state, vars(isF), names(isF), origin(isF));
%                 removed = removed | isF;
%             end
% 
%             % Set up state with remaining variables
%             state = initStateAD@ReservoirModel(model, state, vars(~removed), names(~removed), origin(~removed));
%             % Account for dissolution changing variables
%             if model.disgas
%                 rsSat = model.getProp(state, 'RsMax');
%                 rs = ~st{1}.*rsSat + st{1}.*x;
%                 % rs = rs.*(value(sO) > 0);
%                 state = model.setProp(state, 'rs', rs);
%             end
% 
%             if model.vapoil
%                 rvSat = model.getProp(state, 'RvMax');
%                 rv = ~st{2}.*rvSat + st{2}.*x;
%                 % rv = rv.*(value(sG) > 0);
%                 state = model.setProp(state, 'rv', rv);
%                 % No rv, no so -> zero on diagonal in matrix
%                 bad_oil = value(sO) == 0 & value(rv) == 0;
%                 if any(bad_oil)
%                     sO(bad_oil) = 1 - sW(bad_oil) - value(sG(bad_oil));
%                     state = model.setProp(state, 'sO', sO);
%                 end
%             end
%         end
            
    end
    
end