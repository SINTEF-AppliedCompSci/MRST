function D = computeTOFandTracerAverage(state, G, rock, varargin)
%Executes computeTOFandTracer for a series of states and averages
%
% SYNOPSIS:
%    D  = computeTOFandTracerAverage(state, G, rock, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   Computes time of flight for a series of time steps and avarages.
%
%
% REQUIRED PARAMETERS:
%   G     - Grid structure.
%
%   rock  - Rock data structure.
%           Must contain a valid porosity field, 'rock.poro'.
%
%   state - Reservoir and well solution structure either properly
%           initialized from functions 'initResSol' and 'initWellSol'
%           respectively, or the results from a call to function
%           'solveIncompFlow'.  Must contain valid cell interface fluxes,
%           'state.flux'.
%
%
% OPTIONAL PARAMETERS:
%   wells   - Well structure as defined by function 'addWell'.  May be empty
%             (i.e., wells = []) which is interpreted as a model without any
%             wells.
%   dt      - List of timesteps to use in averaging
%   max_tof - Threshold
%   min_tof - Threshold
%
%   Additional parameters will be passed on to 'computeTOFandTracer', which
%   computes the flow diagnostics. See this function for further documentation.
%
% RETURNS:
%   D - struct that contains the basis for computing flow diagnostics:
%       'inj'     - list of injection wells
%       'prod'    - list of production wells
%       'tof'     - time-of-flight and reverse time-of-flight returned
%                   as an (G.cells.num x 2) vector
%       'itracer' - steady-state tracer distribution for injectors
%       'ipart'   - tracer partition for injectors
%       'ptracer' - steady-state tracer distribution for producers
%       'ppart'   - tracer partition for producers
%       'isubset' - flooding volumes which satisfy the min/max tof
%                   thresholds (if set). Represents a kind of frequency.
%       'psubset' - drainage volumes which satisfy the min/max tof
%                   thresholds (if set). Represents a kind of frequency.
%
% EXAMPLE:
%
% SEE ALSO:
%   `computeTOFandTracer`

%{
Copyright 2009-2018 SINTEF Digital, Applied Mathematics.

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


% Process optional parameters
opt = struct(...
    'wells', [],      ...
    'dt', [],         ...
    'max_tof', [],    ...
    'min_tof', [],    ...
    'diagnostics', [] ...
    );

[opt, extra] = merge_options(opt, varargin{:});

%Find timestep sizes
%if (isfield(state{1}, 'time'))
    %times = cellfun(@(x) x.time, state);
    %times = [0; times];
    %dts = (times(2:end) - times(1:end-1)) / times(end);
if (~isempty(opt.dt))
    dts = opt.dt;
else
    N = numel(state);
    dt = 1/N;
    dts = repmat(dt, N);
end

D = [];

%Accumulate tracer and subsets
h = waitbar(0, 'Step 0');
N = numel(state);
for idx=1:N
    if ~isempty(opt.diagnostics)
        % Use precomptued values
        D_new = opt.diagnostics{idx};
    else
        if (~isempty(opt.wells))
            extra = {'wells', opt.wells{idx}, extra{:}};
        end
        D_new = computeTOFandTracer(state{idx}, G, rock, extra{:});
    end

    if (all(isinf(D_new.tof(:))))
        continue;
    end

    compute_subsets = and(~isempty(opt.min_tof), ~isempty(opt.max_tof));

    %Initialize our average D with zeros
    if (isempty(D))
        D.tof = zeros(size(D_new.tof));
        D.itracer = zeros(size(D_new.itracer));
        D.ptracer = zeros(size(D_new.ptracer));
        D.ipart = zeros(size(D_new.ipart));
        D.ppart = zeros(size(D_new.ppart));

        if (compute_subsets)
            D.isubset = zeros(size(D_new.ipart));
            D.psubset = zeros(size(D_new.ppart));
        end

        num_wells = numel(D_new.inj) + numel(D_new.prod);
        D.inj_avg = zeros(1, num_wells);
        D.prod_avg = zeros(1, num_wells);
    end

    dt = dts(idx);

    %Accumulate TOF and tracer values, etc
    D.tof = D.tof + dt*D_new.tof;
    D.itracer = D.itracer + dt*D_new.itracer;
    D.ptracer = D.ptracer + dt*D_new.ptracer;
    D.inj_avg(D_new.inj) = D.inj_avg(D_new.inj) + dt;
    D.prod_avg(D_new.prod) = D.prod_avg(D_new.prod) + dt;

    if (compute_subsets)
        % Drainage volumes
        data = D.tof(:,2);
        psubset = data >= opt.min_tof & data <= opt.max_tof;

        % Flooding volumes
        data = D.tof(:,1);
        isubset = data >= opt.min_tof & data <= opt.max_tof;

        D.isubset = D.isubset + dt*isubset;
        D.psubset = D.psubset + dt*psubset;
    end

    waitbar(idx/N,h,['Step ', num2str(idx)])
end
delete(h);

%Now, partition the domain
[~,D.ipart] = max(D.itracer,[],2);
[~,D.ppart] = max(D.ptracer,[],2);

%Determine which wells are injectors / producers.
injectors = D.inj_avg > D.prod_avg;
D.inj = find(injectors);
D.prod = find(~injectors);

D.isvalid = true;
end
