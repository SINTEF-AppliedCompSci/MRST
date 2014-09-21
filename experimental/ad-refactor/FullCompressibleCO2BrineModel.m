classdef FullCompressibleCO2BrineModel < ReservoirModel

    % ============================= Class properties =============================
    properties
        
        % Reference temperature at a given reference depth; temperature gradient
        T_ref       % in Kelvin                 
        T_ref_depth % in meters below sea level 
        T_grad      % in degrees Kelvin / 1000 m

        % slope of model
        slope
        slopedir
        
        % custom fluid objects for brine and CO2
        cfluid
        bfluid
        
    end
    % ============================== Public methods ==============================
    methods %(Access = public)

        %% Constructor
        function model = FullCompressibleCO2BrineModel(Gt, rock, tinfo, varargin)
            opt.EOSCO2   = []; 
            opt.EOSBRINE = [];
            opt.verbose  = mrstVerbose();
            opt.slope    = 0;                       % in radians
            opt.slopedir = [1 0];                   % default is slope towards east
            opt.rhoBrine = 1020 * kilogram / meter^3; % if EOS not provided
            opt.mu       = [5.36108e-5, 5.4e-5];    % Default mu values [CO2, brine]
            opt.constantVerticalDensity = false;    % true for semi-compressible model
            opt.nonlinearTolerance = 1e-6;
            
            opt = merge_options(opt, varargin{:});
            
            % Ensuring equations of state are defined
            
            if isempty(opt.EOSCO2)
                opt.EOSCO2 = CO2props('rho_big_trunc','');
            end
            if isempty(opt.EOSBRINE)
                opt.EOSBRINE = FullCompressibleCO2BrineModel.makeConstantEOS(opt.rhoBrine);
            end
            opt.EOSCO2   = FullCompressibleCO2BrineModel.defineAdditionalDerivatives(opt.EOSCO2);
            opt.EOSBRINE = FullCompressibleCO2BrineModel.defineAdditionalDerivatives(opt.EOSBRINE);
            
            % Inherited properties.  The 'fluid' field of ReservoirModel will
            % remain empty, as we use the fields 'cfluid' and 'bfluid' instead.
            model@ReservoirModel(Gt, rock, []);
            model.water = true;
            model.gas   = true;
            model.oil   = false;

            model = model.setupOperators(Gt, rock, varargin{:});
            
            % Other properties
            model.T_ref       = tinfo{1}; 
            model.T_ref_depth = tinfo{2}; 
            model.T_grad      = tinfo{3}; 
            model.slope       = opt.slope;
            model.slopedir    = opt.slopedir;
            model.cfluid = setupFluid(opt.EOSCO2,   opt.mu(1), opt.slope, ...
                                      model.T_grad, opt.constantVerticalDensity);
            model.bfluid = setupFluid(opt.EOSBRINE, opt.mu(2), opt.slope, ...
                                      model.T_grad, opt.constantVerticalDensity);
          
        end
        
        function state = includeComputedCaprockValues(model, state, quick)
        % Given a state with values at the co2-brine interface, compute the
        % corresponding values at the caprock level (pressure and density)
            
        %tempFun  = @(x,i) state.extras.tI(i) - (x - state.h(i) .* cos(model.slope) / 1000;
            cos_theta = cos(model.slope);
            
            if quick
                % 'quick and dirty', using Taylor
                [IEta,~, Eta] = model.cfluid.h_integrals(state.pressure, state.extras.tI);
                
                state.extras.tTop = state.extras.tI - ...
                                    state.h * cos_theta .* model.T_grad / 1000;
                
                pdiff = double(state.h) * norm(gravity) * cos_theta .* double(state.extras.rhoI);
                if ~isempty(IEta)
                    pdiff = pdiff .* double(IEta(-state.h));
                end
                state.extras.pTop   = state.pressure - pdiff;
                state.extras.rhoTop = state.extras.rhoI .* Eta(-state.h);
            else 
                % use integration to get more exact values
                
                num_cells           = numel(double(state.h));
                state.extras.tTop   = zeros(num_cells, 1);
                state.extras.pTop   = zeros(num_cells, 1);
                state.extras.rhoTop = zeros(num_cells, 1);
                
                for i = 1:num_cells
                    tempFun = @(z) state.extras.tI(i) - (z - state.h(i)) * cos_theta / 1000;
                    state.extras.tTop(i) = tempFun(0);
                    state.extras.pTop(i) = numIntTopPress(state.pressure(i), ...
                                                          tempFun, ...
                                                          state.h(i), ...
                                                          model.cfluid.rho, ...
                                                          norm(gravity) * cos_theta);
                end
                state.extras.rhoTop = model.cfluid.rho(state.extras.pTop, state.extras.tTop); 
            end
        
            % @@ If 'quick' has been requested, we could also imagine using
            % Taylor to compute the density at top. 
            state.extras.rhoTop = model.cfluid.rho(state.extras.pTop, state.extras.tTop);
            if ~isempty(IEta)
                
            else
                % Density is either constant in height or always constant
                
            end
        
        end
    end
    % ============================== Private methods ==============================
    methods %(Access = protected)

        function model = setupOperators(model, Gt, rock, varargin)
            
            operators = setupSimCompVe(Gt, rock, varargin{:});
            operators.pv = operators.pv ./ Gt.cells.H; % @@ Necessary since
                                                       % we use a height-formulation.
            model.operators = operators;
        end 
    % ----------------------------------------------------------------------------
        function [problem, state] = ...
                getEquations(model, state0, state, dt, drivingForces, varargin)
            
            [problem, state] = ...
                equationsCO2BrineCompressible(state0, state, dt, model.G, ...
                                              drivingForces, model.operators, ...
                                              model.cfluid, model.bfluid, ...
                                              model.T_ref, model.T_ref_depth, ...
                                              model.T_grad, model.slope, ...
                                              model.slopedir, varargin{:});
        end
    % ----------------------------------------------------------------------------
    function [fn, index] = getVariableField(model, name)
    % Get the index/name mapping for the model (such as where
    % pressure or water saturation is located in state)
       switch(lower(name))
         case {'h', 'height'}
           fn = 'h';
           index = 1;
         % case {'t', 'temperature'}
         %   fn = 'T';
         %   index = 1;
         % case {'sw', 'water'}
         %   index = find(strcmpi(model.saturationVarNames, 'sw'));
         %   fn = 's';
         % case {'so', 'oil'}
         %   index = find(strcmpi(model.saturationVarNames, 'so'));
         %   fn = 's';
         % case {'sg', 'gas'}
         %   index = find(strcmpi(model.saturationVarNames, 'sg'));
         %   fn = 's';
         % case {'s', 'sat', 'saturation'}
         %   index = 1:numel(model.saturationVarNames);
         %   fn = 's';
         case {'pressure', 'p'}
           index = 1;
           fn = 'pressure';
         % case 'wellsol'
         %   % Use colon to get all variables, since the wellsol may
         %   % be empty
         %   index = ':';
         %   fn = 'wellSol';
         otherwise
           % This will throw an error for us
           [fn, index] = getVariableField@PhysicalModel(model, name);
       end
    end
    
        
    % ----------------------------------------------------------------------------    
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok

        %           [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
           
            % computing pressure increment
            dp = dx{problem.indexOfPrimaryVariable('pressure')};
            dp = sign(dp) .* min(abs(dp), abs(model.dpMax .* state.pressure));
            
            % computing height increment
            dh    = dx{problem.indexOfPrimaryVariable('height')};
            dhMax = model.dsMax .* model.G.cells.H;
            dh    = sign(dh) .* min(abs(dh), dhMax);
            
            % computting well-related increments
            dq   = dx{problem.indexOfPrimaryVariable('q')};
            dbhp = dx{problem.indexOfPrimaryVariable('bhp')};
            
            % Updating state with the new incrmeents
            state.pressure     = state.pressure + dp;
            state.h            = state.h + dh;
            if ~isempty(dbhp)
                dbhp = sign(dbhp) .* min(abs(dbhp),  ...
                                         abs(model.dpMax.* vertcat(state.wellSol.bhp)));
                for w = 1:numel(state.wellSol)
                    state.wellSol(w).bhp = state.wellSol(w).bhp + dbhp(w);
                    state.wellSol(w).qGs = state.wellSol(w).qGs + dq;
                    n_ix = state.wellSol.qGs < 0;
                    state.wellSol(w).qGs(n_ix) = 0;
                end
            end
            
            % capping values where necessary
            state.h(state.h<0) = 0;
            ex_ix              = (state.h > model.G.cells.H);
            state.h(ex_ix)     = model.G.cells.H(ex_ix);
            
            % Report is dummy, for now
            report = [];
        end
    end

    % ============================== Static methods ==============================

    methods(Static)
        
        % ----------------------------------------------------------------------------
        function EOS = makeConstantEOS(rho)
            EOS.rho = @(p, t) rho * ones(numel(double(p)), 1);
        end
        
        % ----------------------------------------------------------------------------
       
        function EOS = defineAdditionalDerivatives(EOS)
        % add beta2, gamma2 and chi functions if they are not already there, and if
        % the functions to construct them are avaialble.
            if isfield(EOS, 'rhoDPP') && ~isfield(EOS, 'beta2')
                EOS.beta2  = @(p,t) EOS.rhoDPP(p, t) ./ EOS.rho(p,t);
            end
            if isfield(EOS, 'rhoDTT') && ~isfield(EOS, 'gamma2')
                EOS.gamma2 = @(p,t) EOS.rhoDTT(p, t) ./ EOS.rho(p,t);
            end
            if isfield(EOS, 'rhoDPT') && ~isfield(EOS, 'chi')
                EOS.chi = @(p,t) EOS.rhoDPT(p, t) ./ EOS.rho(p, t);
            end
        end
    end
end

% ====================== Helper functions (not methods) ======================


% ----------------------------------------------------------------------------

function fluid = setupFluid(EOS, mu, slope, Tgrad, vconst)
    fluid.mu  = @(p, t) mu * ones(numel(double(p)), 1); % @@ can be rewritten to allow variable mu
    fluid.rho = @EOS.rho;   % Should be a function of p and t
    
    % define functions to correct for variable vertical density
    fluid.h_integrals = etaIntegralFunctions(EOS, slope, Tgrad, vconst);
end

% ----------------------------------------------------------------------------

function [Ieta, INupEta, Eta] = IetaAndINupEtaAndEta(p, t, EOS, Gct, gct)
    EOS.compressible = 'full'; % required by the etaIntegrals function
    [Ieta, INupEta, ~, ~, Eta] = etaIntegrals(EOS, p , t, Gct, gct); 
end

% ----------------------------------------------------------------------------

function fun = etaIntegralFunctions(EOS, slope, Tgrad, vconst)
    if complete_eos(EOS) && ~vconst
        gct = norm(gravity) * cos(slope);
        Gct = Tgrad / 1000  * cos(slope);
        fun = @(p, t) IetaAndINupEtaAndEta(p, t, EOS, Gct, gct);
    else
        % We do not have a complete, compressible equation of state, or
        % alternatively, the user has requested constant vertical
        % density, so the correction functions are empty
        fun = @(p, t) deal([], [], []);
    end
end

% ----------------------------------------------------------------------------

function res = complete_eos(EOS)
% Check if EOS has all the required functions to be useable for
% approximating vertical density profiles
    
    if isfield(EOS, 'compressible') && ~strcmp(EOS.compressible, 'full')
        res = false;  % the EOS has explicitly been declared not to support
                      % vertical density profiles.  @@ Not particularly
                      % elegant, but a compatibility issue with older code.
        return;
    end

    contains = @(name) isfield(EOS, name) && isa(EOS.(name), 'function_handle');
    
    % Return true if EOS contains all of the following functions:
    res = all(cellfun(contains, {'rho', 'beta', 'gamma', 'chi', 'beta2', 'gamma2'}));
end

% ----------------------------------------------------------------------------

function Ptop = numIntTopPress(Pref, tfun, h, rhofun, g_cos_t)
    res = ode23(@(z, p) g_cos_t * rhofun(p, tfun(z)), [h, 0], Pref);
    Ptop = res.y(end);
end
