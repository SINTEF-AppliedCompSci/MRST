function varargout = runScheduleADIwithVE(initState, G, rock, system, schedule, varargin)
% Given a schedule and a system, solve for all time steps
%
% SYNOPSIS:
%  [wellSols, states, its] = runScheduleADI(initState, G, rock, system, schedule)
%  [wellSols, states, its] = runScheduleADI(initState, G, rock, system, ...
%                                           schedule, 'pn', pv, ...)
% PARAMETERS:
%   initState - The initial state at t = 0;
%
%   G         - A valid grid. See grid_structure.
%
%   rock      - A valid rock structure. Should contain an Nx1 array
%               'poro' containing cell wise porosity values. A permeability
%               field is not *needed* for all the ad-fi solvers as they can
%               work directly with transmissibilities, but it is
%               highly recommended to supply them in either a Nx1 or Nx3
%               array. N is here equal to G.cells.num.
%
%  system     - System configuration as defined by initADISystem.
%
%  schedule   - Schedule (usually found in the deck.SCHEDULE field from
%               the return value of readEclipseDeck from the deckformat
%               module). This fully defines the well configurations for all
%               timesteps.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%           
%   Verbose - If verbose output should be outputted. Defaults to
%             mrstVerbose.
%
%   writeOutput - Save output to the cache folder. This can be practical
%                 when states becomes too big to solve in memory or when
%                 running adjoint simulations.
%
%   outputName  - The string which prefixes .mat files written if
%                 writeOutput is enabled. Defaults to 'state'.
%
%
% RETURNS:
%   wellSols - Well solution struct for each timestep. Cellarray of size
%              Ntx1.
%
%   states (OPTIONAL) - State solution struct for each timestep. Cellarray
%                       of size Ntx1. Note that as this can be come
%                       prohibitively big for long simulations this should
%                       be only outputted if neded.
%
%   its (OPTIONAL) - Nonlinear iteration count for each timestep.
%
%
% SEE ALSO:
%   solvefiADI

%{
#COPYRIGHT#
%}

% $Date: $
% $Revision: $

opt = struct('Verbose', mrstVerbose, ...
             'writeOutput', false,...
             'outputName', 'state');
         
opt = merge_options(opt, varargin{:});

vb = opt.Verbose;
outputStates = nargout > 1;
outputIter = nargout > 2;

%--------------------------------------------------------------------------
dt = schedule.step.val;
tm = cumsum(dt);
dispif(vb, '*****************************************************************\n')
dispif(vb, '********** Starting simulation: %5.0f steps, %5.0f days *********\n', numel(dt), tm(end)/day)
dispif(vb, '*****************************************************************\n')
%--------------------------------------------------------------------------
if opt.writeOutput
    directory = fullfile(fileparts(mfilename('fullpath')), 'cache');
    %delete existing output
    delete(fullfile(directory, [opt.outputName, '*.mat']));
    % output file-names
    outNm  = @(tstep)fullfile(directory, [opt.outputName, sprintf('%05.0f', tstep)]);
end


%--------------------------------------------------------------------------
% output initState
state = initState;
if opt.writeOutput, save(outNm(0), 'state');end;
%--------------------------------------------------------------------------
% collect all wellsols in cell wellsols
wellSols = cell(numel(dt), 1);

if outputStates
    states = cell(numel(dt)+1, 1);
    states{1} = state;
end

iter = zeros(numel(dt),1);

%--------------------------------------------------------------------------
% default is to report all steps
if ~isfield(schedule.step, 'repStep')
    schedule.step.repStep = true(numel(dt), 1);
end

prevControl = -1;
timero = tic;
repStep = 0;

for tstep = 1:numel(dt)
    dispif(vb, 'Time step %5.0f of %d\n', tstep, numel(dt));
    control = schedule.step.control(tstep);
    if control~=prevControl
        if(control==0)
            W=[]; 
        else
            if(~any(strcmp(G.type,'topSurfaceGrid')))
                W = processWells(G, rock, schedule.control(control), 'Verbose', opt.Verbose);
            else
                W3D = processWells(G.parent, rock.parent, schedule.control(control), 'Verbose', opt.Verbose);
                W = convertwellsVE_s(W3D, G.parent, G, rock, 'ip_tpf');
                for i=1:numel(W3D);
                   W(i).bhpLimit=W3D(i).bhpLimit; 
                end
            end
            if isempty(W)
                rock.perm = ones(G.cells.num,1);
                W = addWell([], G, rock, 1, 'Val', 0, 'Type', 'rate', 'sign', 1);
                W.poly = 0;
                W.bhpLimit = 0;
            end
        end
    end

    [state, its] = solvefiADI(state, dt(tstep), W, G, system);
    
    s = [W.sign];
    
    insideLimits = s.*[state.wellSol.bhp] <= [W.bhpLimit].*s;
    if ~all(insideLimits) && opt.Verbose;
        fprintf('Well(s) outside limit: ')
        fprintf('%s ', W(~insideLimits).name);
        fprintf('\n');
    end
    
    iter(tstep) = its;
    wellSols{tstep} = state.wellSol;
    wellSols{tstep} = addWellInfo(wellSols{tstep}, W);
    if outputStates
        states{tstep + 1} = state;
    end
    prevControl = control;
    if opt.writeOutput && schedule.step.repStep(tstep)
        repStep = repStep + 1;
        save(outNm(repStep), 'state'); 
    end
    dispif(~opt.Verbose, 'Step %4g of %4g (Used %3g iterations)\n', tstep, numel(dt), its);
%     dispstr = sprintf('\nStep %4g of %4g\n', tstep, numel(dt));
%     if tstep == 1
%         dispif(~opt.Verbose, dispstr)
%     else
%         dispif(~opt.Verbose, [repmat('\b', 1,numel(dispstr)) dispstr])
%     end
end

dispif(vb, '************Simulation done: %7.2f seconds ********************\n', toc(timero))
varargout{1} = wellSols;
if outputStates
    varargout{2} = states;
end

if outputIter
    varargout{3} = iter;
end
end

function wellSol = addWellInfo(wellSol, W)
%nm = fieldnames(W);
nm = {'name', 'sign'};
for k = 1:numel(nm)
    for wnum = 1:numel(W)
        wellSol(wnum).(nm{k}) = W(wnum).(nm{k});
    end
end
end


