function system = initADISystemVE(input, G, rock, fluid, varargin)
% Generate ADI system. For detailed explanation of the different options,
% see the end of the file.
%
% 
% 
% SYNOPSIS:
%   system = initADISystem(deck, G, rock, fluid)
%   system = initADISystem({'Oil', 'Water'}, G, rock, fluid)
%   system = initADISystem(deck, G, rock, fluid, 'pn', pv,...)
%
% PARAMETERS:
%   input - The first input is either a deck (for autodetection of active
%           phases) or a cell array of strings (for manual configuration of
%           active phases). For instance, if one wants to use fluid data
%           from a case with polymer injection for a reference oil/water
%           simulation.
%
%   G     - Valid grid structure. See grid_structure.
%
%   rock  - rock containing fields perm and poro
%
%   fluid - fluid as defined by initDeckADIFluid;
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   
%
% NOTE
%   
% RETURNS:
%   A system struct
%
% SEE ALSO:
%   

%{
#COPYRIGHT#
%}

opt = struct('tol_mb',              1e-7,...
             'tol_cnv',             1e-3,...
             'tol',                 1e-6,...                 
             'tol_relaxation',      .1, ...
             'maxIts',              25,...
             'Verbose',             mrstVerbose,...
             'IgnoreErrors',        false,...
             'cpr',                 false,...
             'cprEllipticSolver',   [],...
             'cprBlockInvert',      true,...
             'relaxation',          1,...
             'relaxRelTol',         .2,...
             'relaxType',           'dampen',...
             'relaxMax',            .5,...
             'relaxInc',            .1,...
             'pscale',              1,...
             'changeWells',         false,...
             'bhpcontrols',         false,...
             'linesearch',          false,...
             'lineRelTol',          .99,...
             'lineIter',            10,...
             'cprRelTol',           1e-3,...
             'cprType',             'colSum',...
             'podbasis',            [], ...
             'simComponents',       [],...
             'VE',true);

opt = merge_options(opt, varargin{:});

comp = getSupportedComponents();

s = opt.simComponents;
hasDeck = isfield(input, 'RUNSPEC');
if hasDeck
    % We have recieved a deck and should autodetect active phases
    active = @(field) isfield(input.RUNSPEC, upper(field));
    f = fields(comp);
    for i = 1:numel(f)
        comp.(f{i}) = active(f{i});
    end
    if isempty(s)
        s = setupSimComp(G, rock, 'deck', input);
    end
else
    % We have a cell array of components describing the type of system
    f = fields(comp);
    for i = 1:numel(f)
        comp.(f{i}) = any(strcmpi(input, f{i}));
    end
    if isempty(s)
        s = setupSimComp(G, rock);
    end
end


%assert(comp.oil && comp.water, 'Water and Oil phases must be present!')

if isempty(opt.cprEllipticSolver)
    if exist('agmg', 'file') == 2
        % AGMG is a useful elliptic solver. We default to using it if it is
        % installed.
        opt.cprEllipticSolver = @agmg;
    else
        opt.cprEllipticSolver = @mldivide;
    end
end

if comp.gas
    if comp.oil && comp.water,
        if comp.disgas
            dispif(opt.Verbose, 'Found a black oil system...\n')
            system.stepFunction = @(state0, state, meta, dt, W, G, system) stepBlackOil(state0, state, meta, dt, G, W, system, fluid);
            system.getEquations = @eqsfiBlackOilExplicitWells;
            system.cellwise = 1:4;    
            system.cpr.gas = [3 4];
            system.cpr.active = 1:3;
        else
            error('Three phase with not disolved gas not implemented') 
        end
    else
       if(comp.water)
           error('Pure Gas water system not implemented use 3 phase system')
       else
           if(comp.disgas)
               if(opt.VE)
                   dispif(opt.Verbose, 'Found a Oil/Gas with disolution with VE option...\n')
                   system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin)...
                       stepBlackOilOGVE(state0, state, meta, dt, G, W, system, fluid, varargin{:});
                   system.getEquations = @eqsfiBlackOilExplicitWellsOGVE;
                   system.updateFinal  = @(state, state0) updateFinalVE(state, state0, fluid,G);
                   system.cellwise = 1:3;
                   system.cpr.gas = [2 3]
                   system.cpr.active = 1:2;
               else
                   dispif(opt.Verbose, 'Found a Oil/Gas with disolution ...\n')
                   system.stepFunction = @(state0, state, meta, dt, W, G, system) ...
                        stepBlackOilOG(state0, state, meta, dt, G, W, system, fluid);
                   system.getEquations = @eqsfiBlackOilExplicitWellsOG;
                   system.updateFinal  = @(state, state0) updateFinal(state, state0);
                   system.cellwise = 1:3;
                   system.cpr.gas = [2 3]
                   system.cpr.active = 1:2;
                   
               end

               
               
            %error('Oil/Gas system with disolution not implemented');
           else
               if(opt.VE)
                dispif(opt.Verbose, 'Found a Oil/Gas system with VE option...\n')
                system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin)...
                    stepOG(state0, state, meta, dt, G, W, system, fluid, varargin{:});
                system.getEquations = @eqsfiOGExplicitWellsVE;
                system.updateFinal  = @(state, state0) updateFinal(state, state0);
                system.cellwise = 1:2;
                system.cpr.gas = [];
                system.cpr.active = 1:2;  
               else
                 dispif(opt.Verbose, 'Found a Oil/Gas system...\n')
                system.stepFunction = @(state0, state, meta, dt, W, G, system)...
                    stepOG(state0, state, meta, dt, G, W, system, fluid);
                system.getEquations = @eqsfiOGExplicitWells;
                system.updateFinal  = @(state, state0) updateFinal(state, state0);
                system.cellwise = 1:2;
                system.cpr.gas = [];
                system.cpr.active = 1:2;    
               end
           end
       end
    end    
else    
    if comp.polymer
        dispif(opt.Verbose, 'Found a polymer system...\n')
        system.stepFunction = @(state0, state, meta, dt, W, G, system) stepOWPolymer(state0, state, meta, dt, G, W, system, fluid);
        system.getEquations = @eqsfiOWPolymerExplicitWells;
        system.cellwise = 1:3;
        system.cpr.gas = [];
        system.cpr.active = 1:3;
        assert(~comp.gas || opt.IgnoreErrors, 'Gas with polymer not currently supported...')
    else
       dispif(opt.Verbose, 'Found a oil/water system...\n')
        system.stepFunction = @(state0, state, meta, dt, W, G, system) stepOW(state0, state, meta, dt, G, W, system, fluid);
        system.getEquations = @eqsfiOWExplicitWells;
        system.cellwise = 1:2;
        system.cpr.gas = [];
        system.cpr.active = 1:2;
    end
end



% Nonlinear parameters
% Total material balance fpr all cells
system.nonlinear.tolMB          = opt.tol_mb;
% Maximum saturation error
system.nonlinear.tolCNV         = opt.tol_cnv;
system.nonlinear.tol         = opt.tol;

% Maximum number of Newton-Raphson iterations
system.nonlinear.maxIterations  = opt.maxIts;

% Should relaxation be used
system.nonlinear.relaxation     = opt.relaxation;
% The type of relaxation - currently 'dampen' or 'sor' for dampening and
% successive overrelaxation respectively.
system.nonlinear.relaxType      = opt.relaxType;
% The relative tolerance used to determine if a system is oscillating and
% relaxation should be used. A simple test for Newton oscillations is
% performed at each iteration by checking if there are alternating values
% which are approximately equal at iteration i-2 using this tolerance.
system.nonlinear.relaxRelTol    = opt.relaxRelTol;
% Maximum amount of relaxation.
system.nonlinear.relaxMax       = opt.relaxMax;
% Increment of the relaxation parameter each time oscillations are
% detected.
system.nonlinear.relaxInc       = opt.relaxInc;

% Should rate wells be changed to pressure driven when BHP limits are exceeded?
% This option is not well tested.
system.nonlinear.changeWells    = opt.changeWells;
% Force use of BHP controls. If BHP values are defaulted care should be
% taken with this option.
system.nonlinear.bhpcontrols    = opt.bhpcontrols;

% When the Newton iterations have problems with reducing the error, for
% instance near points where the gas saturation is flashed a line search
% can be employed which seeks along the Newton increment vector while
% trying to reduce the error.
system.nonlinear.linesearch  = opt.linesearch;
system.nonlinear.lineRelTol  = opt.lineRelTol;
system.nonlinear.lineIter    = opt.lineIter;

% CPR options
% Use CPR preconditioner if defined for this system
system.nonlinear.cpr            = opt.cpr;
% How the pressure systems are decoupled when preconditioning. Can be
% 'colSum' or 'diag'.
system.nonlinear.cprType        = opt.cprType;
system.nonlinear.cprRelTol      = opt.cprRelTol;
system.nonlinear.cprBlockInvert = opt.cprBlockInvert;

% The specialized solver used for the elliptic pressure-like problem.
system.nonlinear.cprEllipticSolver = opt.cprEllipticSolver;

% Give a list of active components
system.activeComponents = comp;
% Useful constants.
system.s = s;
system.pscale = opt.pscale;

% If an orthogonal basis is provided it can be used for model reduction.
system.podbasis = opt.podbasis;

end


function ss = getSupportedComponents()
   s = { 'oil', 'water', 'gas', 'polymer', 'disgas', 'T', 'MI', 'vapoil' };
    for i = 1:numel(s)
        ss.(s{i}) = false;
    end
end
