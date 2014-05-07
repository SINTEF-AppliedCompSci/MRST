classdef fullCompressibleCO2BrineModel < physicalModel

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
    methods (Access = public)

        %% Constructor
        function model = fullCompressibleCO2BrineModel(Gt, tinfo, varargin)
            opt.EOSCO2   = []; 
            opt.EOSBRINE = [];
            opt.slope    = 0;                    % in radians
            opt.slopedir = [1 0];                % default is slope towards east
            opt.rhoBrine = 1020 kg / meter^3;    % if EOS not provided
            opt.mu       = [5.36108e-5, 6.5e-4]; % Default mu values [CO2, brine]
            op5.constantVerticalDensity = false; % true for semi-compressible model
            opt = merge_options(opt, varargin{:});
            
            % Ensuring equations of state are defined
            if isempty(opt.EOSCO2)
                opt.EOSCO2 = CO2props('rho_big_trunc');
            end
            if isempty(opt.EOSBRINE)
                opt.EOSBRINE = makeConstantEOS(opt.rhoBrine);
            end
            
            % Inherited properties
            model.G         = Gt;
            model.name      = 'fully_compressible_CO2_brine';
            model.operators = model.setupOperators(Gt, rock, varargin{:});
            model.oil       = false;
            model.gas       = true;
            model.water     = true;
            
            % Other properties
            model.T_ref       = tinfo{1}; 
            model.T_ref_depth = tinfo{2}; 
            model.T_grad      = tinfo{3}; 
            model.slope       = opt.slope;
            model.cfluid = setupFluid(opt.EOSCO2,   opt.mu(1), opt.slope, ...
                                      T_grad, opt.constantVerticalDensity);
            model.bfluid = setupFluid(opt.EOSBRINE, opt.mu(2), opt.slope, ...
                                      T_grad, opt.constantVerticalDensity);
          
        end
    end
    % ============================== Private methods ==============================
    methods (Access = private)

        function operators = setupOperators(model, Gt, rock, varargin)
            
            operators = setupSimCompVe(Gt, rock, varargin{:});
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
        function updateState(model, state, problem, dx, drivingForces)
            
        end
    end
    
    % ====================== Helper functions (not methods) ======================
    
    function EOS = makeConstantEOS(rho)
        EOS.rho = @(p, t) rho * ones(numel(double(p)), 1);
    end
    
    % ----------------------------------------------------------------------------
    
    function fluid = setupFluid(EOS, mu, slope, Tgrad, vconst)
        fluid.mu  = @(p, t) mu; % @@ can be rewritten to allow variable mu
        fluid.rho = @EOS.rho;   % Should be a function of p and t
        
        % define functions to correct for variable vertical density
        fluid.h_integrals = etaIntegralFunctions(EOS, slope, Tgrad, vconst);
    end

    % ----------------------------------------------------------------------------
    
    function c = IetaAndINupEta(p, t, EOS, Gct, gct)
        [Ieta, INupEta] = etaIntegrals(EOS, p , t, Gct, gct); 
        c = {Ieta, INupEta};
    end
    
    % ----------------------------------------------------------------------------
    
    function fun = etaIntegralFunctions(EOS, slope, Tgrad, vconst)
        if complete_eos(EOS) && ~vconst
            gct = norm(gravity) * cos(slope);
            Gct = Tgrad / 1000  * cos(slope);
            fun = @(p, t) IetaAndINupEta(p, t, EOS, Gct, gct);
        else
            % We do not have a complete, compressible equation of state, or
            % alternatively, the user has requested constant vertical
            % density, so the correction functions are empty
            fun = @(p, t) {[], []};
        end
    end
    
    % ----------------------------------------------------------------------------
    
    function res = complete_eos(EOS)
    % Check if EOS has all the required functions to be useable for
    % approximating vertical density profiles
        contains = @(name) isfield(EOS, name) && isa(EOS.(name), 'function_handle');
        
        % Return true if EOS contains all of the following functions:
        res = all(contains({'rho', 'beta', 'gamma', 'chi', 'beta2', 'gamma2'}));
    end
    
end
