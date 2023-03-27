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
        
        % @@ hacked
        dpMax
        dsMax
        
        % plot function (if present, called after each ministep)
        plotfun 
        
        % locking pressure for use with indeterminate pressure systems
        locking_pressure
        
    end
    % ============================== Public methods ==============================
    methods %(Access = public)

        %% Constructor
        function model = FullCompressibleCO2BrineModel(Gt, rock, fluid, tinfo, varargin)
            opt.slope    = 0;                       % in radians
            opt.slopedir = [1 0];                   % default is slope towards east
            opt.nonlinearTolerance = 1e-11; % 1e-14
            opt.plotfun = [];
            opt.trans   = [];
            opt.locking_pressure = [];
            
            opt = merge_options(opt, varargin{:});
            
            % Inherited properties.  
            model@ReservoirModel(Gt, rock, fluid, 'nonlinearTolerance', opt.nonlinearTolerance);
            model.water = true;
            model.gas   = true;
            model.oil   = false;

            model = model.setupOperators(Gt, rock, 'trans', opt.trans);

            % Temperature and slope properties
            model.T_ref       = tinfo{1}; 
            model.T_ref_depth = tinfo{2}; 
            model.T_grad      = tinfo{3}; 
            model.slope       = opt.slope;
            model.slopedir    = opt.slopedir;
            
            % Fluid
            model.fluid.gas.h_integrals = model.fluid.gas.h_integrals(model.slope, model.T_grad);
            model.fluid.wat.h_integrals = model.fluid.wat.h_integrals(model.slope, model.T_grad);
            
            % Other stuff
            model.plotfun          = opt.plotfun;
            model.locking_pressure = opt.locking_pressure;

            % @@hack
            model.dpMax = inf;
            model.dsMax = .2;
                       
        end
        
        function state = includeComputedCaprockValues(model, state, quick)
        % Given a state with values at the co2-brine interface, compute the
        % corresponding values at the caprock level (pressure and density)
            
        %tempFun  = @(x,i) state.extras.tI(i) - (x - state.h(i) .* cos(model.slope) / 1000;
           cos_theta = cos(model.slope);
           
           if quick
              % 'quick and dirty', using Taylor
              [IEta, ~, Eta] = model.fluid.gas.h_integrals(state.pressure, state.extras.tI);
              
              state.extras.tTop = state.extras.tI - ...
                  state.h * cos_theta .* model.T_grad / 1000;
              
              pdiff = double(state.h) * norm(gravity) * cos_theta .* double(state.extras.rhoI);
              rho_modif = 1; % default is no modification
              if ~isempty(IEta)
                 pdiff     = pdiff .* double(IEta(-state.h));
                 rho_modif = Eta(-state.h);
              end
              state.extras.pTop   = state.pressure - pdiff;
              state.extras.rhoTop = state.extras.rhoI .* rho_modif;
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
                                                       model.fluid.gas.rho, ...
                                                       norm(gravity) * cos_theta);
              end
              state.extras.rhoTop = model.fluid.gas.rho(state.extras.pTop, state.extras.tTop); 
           end
           
        end
    
        % Following function called by the solver after each ministep.  We can
        % therefore use it to call functions we want to be run after each
        % converged step
        function state = updateAfterConvergence(model, state0, state, dt, drivingForces)
           
        % If user provided a plotting function, call it here
           if ~isempty(model.plotfun)
              model.plotfun(model.G, state, state0, drivingForces);
              drawnow;
              %pause(1);
           end
        end
    end
    
    % ============================== Private methods ==============================
    methods %(Access = protected)

        function model = setupOperators(model, Gt, rock, varargin)
            
           operators = setupSimCompVe(Gt, rock, 'useNewStandard', true, varargin{:});
           operators.pv = operators.pv ./ Gt.cells.H; % @@ Necessary since
                                                      % we use a height-formulation.
           model.operators = operators;
        end 
    % ----------------------------------------------------------------------------
        function [problem, state] = ...
                getEquations(model, state0, state, dt, drivingForces, varargin)
            
            [problem, state] = ...
                equationsCO2BrineCompressible(state0, state, dt, model.G,              ...
                                              drivingForces, model.operators,          ...
                                              model.fluid,                             ...
                                              model.T_ref, model.T_ref_depth,          ...
                                              model.T_grad, model.slope,               ...
                                              model.slopedir,                          ...
                                              'pressure_lock', model.locking_pressure, ...
                                              varargin{:});
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

        %     [state, report] = updateState@ReservoirModel(model, state, problem, dx, drivingForces);
        % @@ Use custom code for the time being, until it is clear how
        % inactive wells should be treated in updateState@ReservoirModel.
           
           % computing pressure increment
            dp = dx{problem.indexOfPrimaryVariable('pressure')};
            dp = sign(dp) .* min(abs(dp), abs(model.dpMax .* state.pressure));
           
            % computing height increment
            dh    = dx{problem.indexOfPrimaryVariable('height')};
            dhMax = model.dsMax .* model.G.cells.H;
            dh    = sign(dh) .* min(abs(dh), dhMax);
           
            % computting well-related increments
            dq   = dx{problem.indexOfPrimaryVariable('qGs')};
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
            pmin = 0;%5.0663e+04;%0;
            state.pressure(state.pressure < pmin) = pmin;
            state.h(state.h < 0) = 0;
            ex_ix                = (state.h > model.G.cells.H);
            state.h(ex_ix)       = model.G.cells.H(ex_ix);
            
            % Report is dummy, for now
            report = [];
        end
    end

    % ============================== Static methods ==============================

    methods(Static)
        
        % ----------------------------------------------------------------------------
       
        % function EOS = add_extras(EOS, slope, Tgrad, vconst)
        % % add beta2, gamma2 and chi functions if they are not already there, and if
        % % the functions to construct them are avaialble.
        %     if isfield(EOS, 'rhoDPP') && ~isfield(EOS, 'beta2')
        %         EOS.beta2  = @(p,t) EOS.rhoDPP(p, t) ./ EOS.rho(p,t);
        %     end
        %     if isfield(EOS, 'rhoDTT') && ~isfield(EOS, 'gamma2')
        %         EOS.gamma2 = @(p,t) EOS.rhoDTT(p, t) ./ EOS.rho(p,t);
        %     end
        %     if isfield(EOS, 'rhoDPT') && ~isfield(EOS, 'chi')
        %         EOS.chi = @(p,t) EOS.rhoDPT(p, t) ./ EOS.rho(p, t);
        %     end
        %     EOS.h_integrals = etaIntegralFunctions(EOS, slope, Tgrad, vconst);            
        % end
    end
end

% ====================== Helper functions (not methods) ======================

function Ptop = numIntTopPress(Pref, tfun, h, rhofun, g_cos_t)
    res = ode23(@(z, p) g_cos_t * rhofun(p, tfun(z)), [h, 0], Pref);
    Ptop = res.y(end);
end

% ----------------------------------------------------------------------------

% function [Ieta, INupEta, Eta] = IetaAndINupEtaAndEta(p, t, EOS, Gct, gct)
%     EOS.compressible = 'full'; % required by the etaIntegrals function
%     [Ieta, INupEta, ~, ~, Eta] = etaIntegrals(EOS, p , t, Gct, gct); 
% end

% ----------------------------------------------------------------------------

% function fun = etaIntegralFunctions(EOS, slope, Tgrad, vconst)
%     if complete_eos(EOS) && ~vconst
%         gct = norm(gravity) * cos(slope);
%         Gct = Tgrad / 1000  * cos(slope);
%         fun = @(p, t) IetaAndINupEtaAndEta(p, t, EOS, Gct, gct);
%     else
%         % We do not have a complete, compressible equation of state, or
%         % alternatively, the user has requested constant vertical
%         % density, so the correction functions are empty
%         fun = @(p, t) deal([], [], []);
%     end
% end

% ----------------------------------------------------------------------------

% function res = complete_eos(EOS)
% % Check if EOS has all the required functions to be useable for
% % approximating vertical density profiles
    
%     if isfield(EOS, 'compressible') && ~strcmp(EOS.compressible, 'full')
%         res = false;  % the EOS has explicitly been declared not to support
%                       % vertical density profiles.  @@ Not particularly
%                       % elegant, but a compatibility issue with older code.
%         return;
%     end

%     contains = @(name) isfield(EOS, name) && isa(EOS.(name), 'function_handle');
    
%     % Return true if EOS contains all of the following functions:
%     res = all(cellfun(contains, {'rho', 'beta', 'gamma', 'chi', 'beta2', 'gamma2'}));
% end

% ----------------------------------------------------------------------------

