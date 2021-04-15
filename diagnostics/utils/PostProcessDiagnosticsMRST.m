function [d] = PostProcessDiagnosticsMRST(problem,varargin)
%Launch flow diagnostics postprocessor GUI for MRST simulation output.
%
% SYNOPSIS:
%   [d] = PostProcessDiagnosticsMRST(problem,varargin)
% 
% DESCRIPTION:
%   Takes the output of an MRST simulation problem and launches the flow
%   diagnostics postprocessor GUI. The simulation should already
%   have been run as PostProcessDiagnosticsMRST.m will only look for 
%   exisiting output.
%
% REQUIRED PARAMETERS:
%  problem - An MRST simulation problem defined using
%            packSimulationProblem.m.
%
% OPTIONAL PARAMETERS:
%   steps - Array of timesteps to be displayed in the GUI.
%
%   maxTOF - Maximum TOF value. Default value is 500 years.
%
%   cleanup - true or false. Delete all existing flow diagnostics 
%             results from previous runs. Default is false.
% 
%   precompute - true or false. Precompute flow diagnostics for all
%                timesteps. Default is true.
%
%   startdate - optionally specify a start date of the format 
%                [day month year] e.g. [1 1 2000] for the 1st January
%                2000.
%
%
% RETURNS:
%   d      - Handle to PostProcessDiagnostics object.
%
% SEE ALSO:
%   `PostProcessDiagnosticsECLIPSE`.

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

opt = struct('style', 'default', ...
    'steps',        [], ...
    'maxTOF', 500*year, ...
    'cleanup',   false, ...
    'precompute', true, ...
    'startdate',    [0 0 0]);

mrstModule add mrst-gui diagnostics deckformat ad-props coarsegrid
mrstVerbose(true)

opt = merge_options(opt, varargin{:});
pth = problem.OutputHandlers.states.getDataPath();
wsHandler = problem.OutputHandlers.wellSols;


% precompute options
precompDir = fullfile(pth, 'mrst_diagnostics');
if opt.cleanup
    cleanupDialogue(precompDir);
end


if opt.precompute
    dd = dir(precompDir);
    %if exist(precompDir,'dir')~=7
    %    mkdir(precompDir);
    %end
    % only precompute if directory is empty or non-existent
    if ~any(~strncmp({dd.name}, '.', 1))
        precomputeDialogue(problem, precompDir);
    end
end

info.date = opt.startdate;
info.time = cumsum(problem.SimulatorSetup.schedule.step.val(:))./day;

% Select which time-steps to include
if ~isempty(opt.steps)
    steps = opt.steps;
else
    steps = uiPreSelectTimeSteps(info);
end

% Setup data for selected steps and load/compute diagnostics
d.maxTOF = opt.maxTOF;
if opt.precompute
    precomp = getPrecomputedDiagnostics(problem, steps);
else
    precomp = [];
end
[d.G, d.Data, d.Gs, valid_ix] = readAndPrepareForPostProcessorMRST(problem, steps, info, precomp);

if ~any(valid_ix)
    return;
end
if ~all(valid_ix) && ~isempty(precomp)
    precomp = precomp(valid_ix);
end

d.Data.wsdata = readwellSolDataForPostProcessor(problem, 'wellSolFields',{'bhp','qWs','qGs','qOs'}', ...
    'startdate', opt.startdate);

d.Data.ws = getPackedSimulatorOutput(problem);

d = PostProcessDiagnostics(d,precomp,'style',opt.style);

