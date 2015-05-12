function D = computeTOFandTracerAverage(state, G, rock, varargin)
%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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
    'wells', [],...
    'dt', [],...
    'max_tof', [], ...
    'min_tof', [] ...
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
    N = numel(W);
    dt = 1/N;
    dts = repmat(dt, N);
end

D = [];

%Accumulate tracer and subsets
h = waitbar(0, 'Step 0');
N = numel(state);
for idx=1:N
    if (~isempty(opt.wells))
        extra = {'wells', opt.wells{idx}, extra{:}};
    end
    D_new = computeTOFandTracer(state{idx}, G, rock, extra{:});
    
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
        D.inj_avg = zeros(num_wells, 1);
        D.prod_avg = zeros(num_wells, 1);
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