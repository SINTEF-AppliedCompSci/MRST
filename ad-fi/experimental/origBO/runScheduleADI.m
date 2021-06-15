function varargout = runScheduleADI(initState, G, rock, system, schedule, varargin)
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


default_outputDir = fullfile(fileparts(mfilename('fullpath')), 'cache');

opt = struct('Verbose', mrstVerbose, 'writeOutput', false, 'outputName', 'state', 'scaling', [], ...
             'startAt', 1, 'outputDir', default_outputDir, 'Wext',[]);

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
    outNm  = @(tstep)fullfile(opt.outputDir, [opt.outputName, sprintf('%05.0f', tstep)]);
end


%--------------------------------------------------------------------------
% output initState
if opt.startAt == 1
    state = initState;
    if opt.writeOutput, save(outNm(0), 'state');end
else
    rf = load(outNm(opt.startAt-1));
    state = rf.state;
end

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
W = [];

for tstep = opt.startAt:numel(dt)
    dispif(vb, 'Time step %5.0f of %d\n', tstep, numel(dt));
    control = schedule.step.control(tstep);
    if control~=prevControl
        if(control==0)
            W=[];
        else
            if(~any(strcmp(G.type,'topSurfaceGrid')))
                W = processWells(G, rock, schedule.control(control), 'Verbose', opt.Verbose);
                if(~isempty(opt.Wext))
                    for i=1:numel(W)
                     if(isfield(opt.Wext,'T'))
                        W(i).T=opt.Wext(i).T;
                     end
                     if(isfield(opt.Wext,'I'))
                        W(i).I=opt.Wext(i).I;%ones(1,size(state.I,2));
                     end
                    end
                end
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
            for k = 1:numel(W)
                if ~isfield(W(k), 'bhpLimit')
                    if W(k).sign == -1
                        W(k).bhpLimit = 0.;
                    elseif W(k).sign == 1
                        W(k).bhpLimit = 1.0e9; % dummy value!!
                    end
                end
            end
        end
    end

    not_converged = true;
    dispif(vb, sprintf('Time step length: %g day.\n', convertTo(dt(tstep), day)))
    state0 = state;

    [state, its] = solvefiADI(state, dt(tstep), W, G, system);

    if(~isempty(W))
    s = [W.sign];

    insideLimits = s.*[state.wellSol.bhp] <= [W.bhpLimit].*s;
    if ~all(insideLimits) && opt.Verbose;
        fprintf('Well(s) outside limit: ')
        fprintf('%s ', W(~insideLimits).name);
        fprintf('\n');
    end
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
        save(outNm(tstep), 'state', 'W');
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


