function [controls] = initControls(schedule, varargin)
% initControls -- Initialize control structure based on well schedule
%
% SYNOPSIS:
%   controls = initControls(schedule, 'pn', pv, ...)
%
% DESCRIPTION:
%   Initialize controls
%
% PARAMETERS:
%   schedule       - schedule structure
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - ControllableWells  :  indices to wells which will be
%                                       controlled (default all)
%               - BHPMaxMin          :  max/min value for bhp-wells,(default [-Inf Inf])
%               - RateMaxMin         :  max/min value for rate-wells,(default [-Inf Inf])
%               - MaxMin             :  matrix of size numControls x 2,
%                                       where each row contains the max and min
%                                       value for the corresponding control. Aternative
%                                       to the two above
%               - NumControlSteps    :  number of control steps (default number of time steps)
%               - LinEqConst         :  linear equality constraints of the
%                                       form Au = b given in the form {A_1,
%                                       b_1, A_2, b_2, ...}
%               - Verbose            :  Display Control vars using dispControls
%
% RETURNS:
%   controls   - Initialized control structure having fields
%                   - well  : #controllable wells x 1 structure having
%                             fields:
%                       - wellNum    : index to current well in schedule (and W)
%                       - values     : values for each timeStep
%                       - bhpMaxMin  :
%                       - rateMAxMin :
%                   - linEqConst : linear equality contraints structure having fields:
%                       - A
%                       - b
%
% SEE ALSO:
%              `initSchedule`

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


numWells    = numel( schedule(1).names );
opt = struct('ControllableWells',   (1:numWells)', ...
             'BHPMinMax',           [-Inf Inf], ...
             'RateMinMax',          [-Inf Inf], ...
             'MinMax',              [], ...
             'NumControlSteps',     numel(schedule), ...
             'LinIneqConst',        [], ...
             'LinEqConst',          [], ...
             'Verbose',             mrstVerbose );

opt = merge_options(opt, varargin{:});
verbose = opt.Verbose;
cw      = opt.ControllableWells;
bhpmm   = opt.BHPMinMax;
ratemm  = opt.RateMinMax;
minMax  = opt.MinMax;
linec   = opt.LinIneqConst;
lec     = opt.LinEqConst;

%--------------------------------------------------------------------------

numControlWells  = length(cw);
wellTypes        = vertcat( schedule(1).types );
controlTypes     = wellTypes(cw);

if isempty( minMax ), minMax = ones(numControlWells, 1)*[-Inf Inf]; end
vals    = [schedule(:).values]';

for k = 1:length(cw)
    well(k).wellNum     = cw(k);
    well(k).values      = vals(:, cw(k));

    type                = controlTypes{k};
    well(k).type        = type;
    if strcmp(type, 'rate')
         well(k).minMax  = [ max( minMax(k,1), ratemm(1) ), ...
                             min( minMax(k,2), ratemm(2) )];
    elseif strcmp(type, 'bhp')
        well(k).minMax  = [ max( minMax(k,1), bhpmm(1) ), ...
                            min( minMax(k,2), bhpmm(2) )];
    end
end

if ~isempty(lec)
    numEq = numel(lec)/2;
    A = []; b = [];
    for k = 1: numEq
        A = [A; lec{2*k-1}];
        b = [b; lec{2*k}];
    end
    controls.linEqConst.A  = A;
    controls.linEqConst.b  = b;
else
    controls.linEqConst = [];
end

if ~isempty(linec)
    numEq = numel(lec)/2;
    %A = []; b = [];
    %for k = 1: numEq
    %    A = [A; lec{2*k-1}];
    %    b = [b; lec{2*k}];
    %end
    %controls.linEqConst.A  = A;
    %controls.linEqConst.b  = b;
else
    controls.linIneqConst = [];
end


controls.well          = well;

if verbose
    dispControls(controls, schedule);
end

controls.numControlSteps = opt.NumControlSteps;
