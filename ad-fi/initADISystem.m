function system = initADISystem(input, G, rock, fluid, varargin)
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
%   fluid - fluid as defined by initDeckADIFluid
%
% OPTIONAL PARAMETERS:
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


opt = struct('tol_mb',              1e-7,...
             'tol_cnv',             1e-3,...
             'tol',                 1e-4,...
             'tol_relaxation',      .1, ...
             'use_ecltol',          true,...
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
             'allowControlSwitching', true, ...
             'allowWellSignChange', false, ...
             'allowCrossFlow', false, ...
             'linesearch',          false,...
             'lineRelTol',          .99,...
             'lineIter',            10,...
             'cprRelTol',           1e-3, ...
             'agmgTol',             1e-10, ...
             'cprType',             'colSum',...
             'podbasis',            [], ...
             'simComponents',       []);

opt = merge_options(opt, varargin{:});


fn = reshape(getSupportedComponents(), [], 1);
s  = opt.simComponents;

if isstruct(input) && isfield(input, 'RUNSPEC') && isstruct(input.RUNSPEC),

   active = isfield(input.RUNSPEC, upper(fn));

   if isempty(s), s = setupSimComp(G, rock, 'deck', input); end

else
   assert (iscellstr(input), ...
          ['First parameter ''input'' must be deck structure or ', ...
           'cell array of strings.']);

   % We have a cell array of components describing the type of system

   input  = reshape(input, [], 1);
   [i, j] = blockDiagIndex(numel(input), numel(fn));
   active = accumarray(j, strcmpi(input(i), fn(j))) > 0;

   if isempty(s), s = setupSimComp(G, rock); end
end

comp = cell2struct(num2cell(active), fn, 1);

%assert(comp.oil && comp.water, 'Water and Oil phases must be present!')

if isempty(opt.cprEllipticSolver)
    if exist('agmg', 'file') == 2
        % AGMG is a useful elliptic solver. We default to using it if it is
        % installed.
        opt.cprEllipticSolver = @(A, B)(agmg(A, B, [], 1e-10));
    else
        opt.cprEllipticSolver = @mldivide;
    end
end
fluid = assignRelPerm(fluid);
if comp.gas
    if comp.oil && comp.water,
        if comp.disgas, ld = 'live'; else ld = 'dead';end
        if comp.vapoil, wd = 'wet' ; else wd = 'dry' ;end
        dispif(opt.Verbose, 'Found a three-phase system: %s oil, %s gas ...\n', ld, wd)
        system.stepFunction = @(state0, state, meta, dt, W, G, system,varargin) stepVO(state0, state, meta, dt, G, W, system, fluid);
        system.getEquations = @eqsfiVO;
        system.cellwise   = 1:3;
        system.cpr.active = 1:3;
        system.cpr.gas    = 3;
%         if comp.disgas
%             dispif(opt.Verbose, 'Found a black-oil system...\n')
%             system.stepFunction = @(state0, state, meta, dt, W, G, system) stepBlackOil(state0, state, meta, dt, G, W, system, fluid);
%             system.getEquations = @eqsfiBlackOil;
%             system.cellwise = 1:4;
%             system.cpr.gas = [3 4];
%             system.cpr.active = 1:3;
%         else
%             error('Three phase with not disolved gas not implemented')
%         end
    else
       if(comp.water)
           error('Pure Gas water system not implemented use 3 phase system')
       else
           if(comp.disgas)
               dispif(opt.Verbose, 'Found a Oil/Gas with disolution...\n')
            system.stepFunction = @(state0, state, meta, dt, W, G, system, varagin) stepBlackOilOG(state0, state, meta, dt, G, W, system, fluid);
            system.getEquations = @eqsfiBlackOilExplicitWellsOG;
            system.updateFinal  = @(state, state0) updateFinal(state, state0);
            system.cellwise = 1:3;
            system.cpr.gas = [2 3];
            system.cpr.active = 1:2;



            %error('Oil/Gas system with disolution not implemented');
           else
               if(~comp.T)

                    dispif(opt.Verbose, 'Found a Oil/Gas system...\n')
                    system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin) stepOG(state0, state, meta, dt, G, W, system, fluid);
                    system.getEquations = @eqsfiOGExplicitWells;
                    system.updateFinal  = @(state, state0) updateFinal(state, state0);
                    system.cellwise = 1:2;
                    system.cpr.gas = [];
                    system.cpr.active = 1:2;
               else
                  if(~comp.MI)
                    dispif(opt.Verbose, 'Found a Oil/Gas system with tempratrue...\n')
                    system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin) stepOGT(state0, state, meta, dt, G, W, system, fluid);
                    system.getEquations = @eqsfiOGTExplicitWells;
                    system.updateFinal  = @(state, state0) updateFinal(state, state0);
                    system.cellwise = 1:3;
                    system.cpr.gas = [];
                    system.cpr.active = 1:3;
                  else
                    dispif(opt.Verbose, 'Found a Oil/Gas system with tempratrue...\n')
                    system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin) stepOGTMI(state0, state, meta, dt, G, W, system, fluid);
                    system.getEquations = @eqsfiOGTMIExplicitWells;
                    system.updateFinal  = @(state, state0) updateFinal(state, state0);
                    system.cellwise = 1:3;
                    system.cpr.gas = [];
                    system.cpr.active = 1:3;
                  end
               end
           end
       end
    end
else
    if comp.polymer,
        dispif(opt.Verbose, 'Found a polymer system...\n');

        system.stepFunction = @(state0, state, meta, dt, W, G, system, varargin) ...
           stepOWPolymer(state0, state, meta, dt, G, W, system, fluid);

        system.getEquations = @eqsfiOWPolymerExplicitWells;
        system.cellwise     = 1:3;
        system.cpr.gas      = [];
        system.cpr.active   = 1:3;

        assert(~comp.gas || opt.IgnoreErrors, ...
               'Gas with polymer not currently supported...');
    else
        [use_T, use_MI] = deal(comp.T, comp.MI);

        dispif(          opt.Verbose, 'Found an oil/water system');
        dispif(use_T  && opt.Verbose, ' with temperature');
        dispif(use_MI && opt.Verbose, ' with minerals');
        dispif(          opt.Verbose, ' ...\n');

        system.stepFunction = @(state0, state, meta, dt, W, G, system,varargin) ...
            stepOW(state0, state, meta, dt, G, W, system, fluid, ...
            'temperature', use_T, 'minerals', use_MI);
        if use_T || use_MI
            system.cellwise     = [1:2,6];
            system.cpr.active   = [1:2,6];
            system.getEquations = @eqsfiOWExplicitWells;
        else
            system.cellwise     = 1:2;
            system.cpr.active   = 1:2;
            system.getEquations = @eqsfiOW;
        end
        
        system.cpr.gas      = [];
        
    end
end



% Nonlinear parameters
% Total material balance fpr all cells
system.nonlinear.tolMB          = opt.tol_mb;
% Maximum saturation error
system.nonlinear.tolCNV         = opt.tol_cnv;
% Tolerance used for simple convergence
system.nonlinear.tol            = opt.tol;
% Use eclipse type tolerance
system.nonlinear.use_ecltol               = opt.use_ecltol;
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

system.well.allowControlSwitching = opt.allowControlSwitching;
system.well.allowWellSignChange = opt.allowWellSignChange;
system.well.allowCrossFlow = opt.allowCrossFlow;

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

%Newton step options
% cap pressure changes (relative)
system.stepOptions.dpMax  = inf;
% cap saturation changes (absolute)
system.stepOptions.dsMax  = .2;
% cap rs changes (relative)
system.stepOptions.drsMax = inf;
% solve well eqs based on updated res props rather than linearized update
system.stepOptions.solveWellEqs = false;
% Number of iterations in solve well eqs based on updated res props
system.stepOptions.maxitWell = 20;
% Number of iterations in solving local equations
system.stepOptions.maxitLocal = 20;


% If an orthogonal basis is provided it can be used for model reduction.
system.podbasis = opt.podbasis;

end

%--------------------------------------------------------------------------
function s = getSupportedComponents()
   s = { 'oil', 'water', 'gas', 'polymer', 'disgas', 'vapoil', 'T', 'MI' };
end
%function s = getSupportedComponents()
%   s = { 'oil', 'water', 'gas', 'polymer', 'disgas', 'T', 'MI' };
%end
