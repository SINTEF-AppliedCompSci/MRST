function [schedule] = initSchedule(W, varargin)
% initSchedule -- Initialize schedule structure based on well W.
%
% SYNOPSIS:
%   shedule = initSchedule(W, 'pn', pv, ...)
%
% DESCRIPTION:
%   Initialize schedule
%
% PARAMETERS:
%   W       - well structure
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - NumSteps  :  number of simulation time steps (default 1)
%               - TotalTime :  total simualtion time (default 1)
%               - TimeSteps :  endtime for each time step assuming t_0 = 0 (alterative to two previous pns)
%               - Verbose   :  display schedule with dispSchedule
%
% RETURNS:
%   schedule   - Initialized numSteps x 1 rate schedule structure having fields
%                   - timeInterval     -- [startTime endTime]
%                   - names            -- {name_1, ..., name_n}
%                   - types            -- {welltype_1, ... , welltype_n}
%                   - values           -- {val_1, ..., val_2}
%
%
% SEE ALSO:
%              `dispSchedule`

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


opt = struct('Verbose',   mrstVerbose, ...
             'NumSteps',  1, ...
             'TotalTime', 1, ...
             'TimeSteps', []);
opt = merge_options(opt, varargin{:});

verbose      = opt.Verbose;
timeSteps    = opt.TimeSteps;

if isempty(timeSteps)
    if ~isfield(opt, 'NumSteps') && ~isfield(opt, 'TotalTime')
        error('Non compatible input ...')
    else
        numSteps = opt.NumSteps;
        totTime  = opt.TotalTime;
        timeSteps = (totTime/numSteps)*(1 : numSteps);
    end
end

ts = timeSteps(:);
intervals = [ [0; ts(1:end-1)] ts(1:end) ];

for k = 1 : length(timeSteps)
    schedule(k).timeInterval    = intervals(k, :);
    schedule(k).names           = {W(:).name}';
    schedule(k).types           = {W(:).type}';
    schedule(k).values          = [W(:).val]';
end

if verbose, dispSchedule(schedule); end
