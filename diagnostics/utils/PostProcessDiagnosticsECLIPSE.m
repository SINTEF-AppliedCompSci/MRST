function [d] = PostProcessDiagnosticsECLIPSE(varargin)
%Launch flow diagnostics postprocessor GUI for ECLIPSE simulation output.
%
% SYNOPSIS:
%   [d] = PostProcessDiagnosticsECLIPSE(varargin)
%
% DESCRIPTION:
%   Takes the output of an ECLIPSE simulation  and launches the flow
%   diagnostics postprocessor GUI. Optionally the path to the .EGRID file
%   can be passed as an input parameter. If no parameters are given a file
%   dialogue box will be launched to select the .EGRID file.
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
% RETURNS:
%   d      - Handle to PostProcessDiagnostics object.
%
% SEE ALSO:
%   `PostProcessDiagnosticsMRST`.

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
    'precompute', true);
mrstModule add mrst-gui diagnostics deckformat ad-props coarsegrid
mrstVerbose(true)

if mod(numel(varargin), 2) == 1  % file-name provided
    filenm = varargin{1};
    [pth, nm, ext] = fileparts(filenm);
    if ~strcmp(ext, 'EGRID')
        filenm = fullfile(pth, [nm, '.', 'EGRID']);
    end
    [opt, extraOpt] = merge_options(opt, varargin{2:end});
else
    [fn, pth] = uigetfile('*.EGRID', 'Select ECLIPSE output grid file (EGRID)');
    if ~isempty(fn) && ischar(fn)
        filenm = fullfile(pth, fn);
    else
        return;
    end
    [opt, extraOpt] = merge_options(opt, varargin{:});
end
assert(exist(filenm, 'file')>0, sprintf('Unable to find file %s', filenm));
[pth, nm] = fileparts(filenm);
casenm = fullfile(pth, nm);

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
        precomputeDialogue(casenm, precompDir);
    end
end
% Select which time-steps to include
rsspec = processEclipseRestartSpec(casenm, 'all');
if ~isempty(opt.steps)
    steps = opt.steps;
else
    steps = uiPreSelectTimeSteps(rsspec);
end

% Setup data for selected steps and load/compute diagnostics
d.maxTOF = opt.maxTOF;
if opt.precompute
    precomp = getPrecomputedDiagnostics(casenm, steps);
else
    precomp = [];
end
[d.G, d.Data, d.Gs, valid_ix] = readAndPrepareForPostProcessor(casenm, steps, rsspec, precomp);
if ~any(valid_ix)
    return;
end
if ~all(valid_ix) && ~isempty(precomp)
    precomp = precomp(valid_ix);
end

try
    d.Data.summary = readEclipseSummaryUnFmt(casenm);
catch
    d.Data.summary = [];
    warning('Not able to read summary for selected case ...\n')
end

d = PostProcessDiagnostics(d,precomp,'style',opt.style, extraOpt{:});

