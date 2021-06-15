function phi = basisOW(states, varargin)%, efrac_p, efrac_s, n_p, n_s)
    % efrac_(s/p) fraction of energy to be included (upper bound) for
    % saturation or pressure
    % n_(s/p) upper bound for the number of eigenvectors to be included

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

    opt = struct('energyfraction_pressure', 1,...
                 'energyfraction_saturation', 1,...
                 'maxno_pressure',          inf,...
                 'maxno_saturation',        inf,...
                 'cluster_n_pressure',      0, ...
                 'cluster_n_saturation',    0, ...
                 'cluster_iter',            1);

    opt = merge_options(opt, varargin{:});



    efrac_p = opt.energyfraction_pressure;
    efrac_s = opt.energyfraction_saturation;
    n_p = opt.maxno_pressure;
    n_s = opt.maxno_saturation;

    nw = numel(states{end}.wellSol);
    nc = numel(states{end}.pressure);

    % Ignore initial state as it has no well solution
    states = horzcat(states{2:end});

    phi.basis = cell(5,1);
    phi.bar = cell(5,1);

    % Pressure
    if isnan(n_p)
        disp 'Including pressure equations in full!'
        basis = speye(nc);
        xbar   = zeros(nc,1);
    else
        pbasis = horzcat(states(:).pressure);
        if opt.cluster_n_pressure
            disp 'Clustering pressure...'
            pbasis = double(yael_kmeans(single(pbasis), opt.cluster_n_pressure, ...
                    'redo', opt.cluster_iter));
        end
        [basis xbar] = createPod(pbasis, efrac_p, n_p);
    end

    phi.basis{1} = basis;
    phi.bar{1} = xbar;


    totsat = horzcat(states(:).s);
    % extract first saturation vector
    if isnan(n_s)
        disp 'Including saturation equations in full!'
        basis = speye(nc);
        xbar   = zeros(nc,1);
    else
        sbasis = totsat(:, 1:2:end-1);
        if opt.cluster_n_saturation
            disp 'Clustering saturation...'
            % Requires some kmeans implementation. Current setup uses yael
            % kmeans which can be found at https://gforge.inria.fr/
            sbasis = double(yael_kmeans(single(sbasis), opt.cluster_n_saturation, ...
                    'redo', opt.cluster_iter));
        end
        [basis xbar] = createPod(sbasis, efrac_s, n_s);
    end

    phi.basis{2} = basis;
    phi.bar{2} = xbar;

    % 1: pressure, oil
    % 2: saturation water
    % rest wells and bc

    % qWs
    for i = 3:5
        % Always include well equations in full
        phi.basis{i} = eye(nw);
        phi.bar{i}   = zeros(nw,1);
    end

%     phi = horzcat(phi{:});
end
