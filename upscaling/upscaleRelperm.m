function [sat_mat, kr, perm, krK] = upscaleRelperm(G, rock, fluid, dp_scale, sat_vec, varargin)
% Upscale relative permeability.
%
% SYNOPSIS:
%   [sat_mat, kr, perm, krK] = upscaleRelperm(G, rock, fluid, dp_scale, sat_vec)
%   [sat_mat, kr, perm, krK] = upscaleRelperm(G, rock, fluid, dp_scale, 'pn', pv, ...)
%
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure, or
%             alternatively a periodic grid as created by
%             makePeriodicGridMulti3d. If using a periodic grid the
%             keyword 'periodic' must be set to true
%
%   rock    - Valid rock data structure.
%
%   fluid   - Fluid object as defined by initSWOFFluidJfunc.
%
%   dp_scale - Pressure drop used for upscaling boundary conditions.
%
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%           - psolver -- function handle to pressure sovler which uses the
%                        interface @(state, Grid, Fluid, BCP, BC, Rock).
%
%           - periodic -- Use periodic grid for upscaling.
%
%           - dir -- Upscale direction
%
%           - dp_vec -- Vector of pressure differentials for boundary
%                       conditions.
%
% RETURNS:
%   sat_mat - Matrix of saturation values. One row for each pressure
%             differential.
%
%   kr      - Relperm values. One for each component in the fluid.
%
%   perm    - Permeability upscaled using pure fluid.
%
%   krK      - Effective permeability. The upscaling is based on first
%              finding an upscaled permeability with a fluid with unitary
%              values (fluid_pure). Once this is found, relperm is found by
%              solving the expression kr * K = krK. Depending on the
%              permeability type (full tensor / diagonal) this is be done
%              by inversion of K or division.
%
% COMMENTS:
%
%
% SEE ALSO:
%

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


opt = struct('psolver' ,   [], ...
             'periodic',   false, ...
             'dir'    ,    1, ...
             'bc_faces',   [],...
             'dp_vec',     []...
             );
opt         = merge_options(opt, varargin{:});
d           = G.griddim;


if isempty(opt.dp_vec)
    opt.dp_vec = dp_scale;
end

psolver = opt.psolver;
if isempty(psolver)
   IP = @(Grid, Rock) computeMimeticIPGp(G, Grid, Rock, opt.periodic);
   psolver = @(state, Grid, Fluid, BCP, BC, Rock) ...
      incompMimetic(state, Grid, IP(Grid, Rock), ...
                    Fluid, 'bcp', BCP, 'bc', BC);
end

periodic    = opt.periodic;
L           = max(G.faces.centroids)-min(G.faces.centroids);
fluid_pure  = initSingleFluid('mu' ,1, 'rho', 1);
dir         = opt.dir;
pv          = poreVolume(G,rock);
N_sat       = numel(sat_vec);
N_dp        = numel(opt.dp_vec);
N_ph        = 2;

% Pairs of opposing boundary condition keywords
bcsides = {'XMin', 'XMax';...
           'YMin', 'YMax';...
           'ZMin', 'ZMax'};


if periodic
    bcl   = cell(d,1);
    bcr   = cell(d,1);
    dp    = cell(d,1);

    for j = 1:d
        bcl{j} = pside([], G, bcsides{j, 1}, 0);
        bcr{j} = pside([], G, bcsides{j, 2}, 0);
        dp{j} = 0;
    end

    dp{1} = dp_scale;
    % Create a periodic grid which couples the edges so that bcl{i}.faces
    % are coupled with bcr{i}.faces for i \in G.griddim

    [Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dp);

    % Define a pressure solver based on the earlier function which uses
    % periodic boundary conditions
    psolver = @(State, Grid, Fluid, BCP, Rock) psolver(State, Grid, Fluid, BCP, [], Rock);

    % Pressure solver for the steady state
    psolver_steady = @(State, Fluid, BCP) psolver(State, Gp, Fluid, BCP, rock);

    % Transmissibilities
    Trans = computeTransGp(G, Gp, rock);

    % Steady state solver
    steadystate = @(state, dt, dummy) simulateToSteadyState(state, Gp, rock, fluid, dt, 'bcp', bcp, 'psolver', psolver_steady, 'trans', Trans);
    upscale = @(Rock, BCP) upscalePermeabilityPeriodic(Gp, BCP, dp_scale, psolver, fluid_pure, Rock, L);

    clear bcl bcr
else
    if isempty(opt.bc_faces);
        % Not supplied, using xMin and xMax
        t1 = pside([],G,'XMin',0);
        t2 = pside([],G,'XMax',0);
        opt.bc_faces{1} = t1.face;
        opt.bc_faces{2} = t2.face;
        clear t1 t2
    end
    assert(numel(opt.bc_faces)==2);

    % Fixed grid, will not create a periodic grid.
    Gp = G;

    % Define pressure solver which uses normal BC
    psolver        = @(State, Grid, Fluid, BC, Rock) psolver(State, Grid, Fluid, [], BC, Rock);

    % Pressure solver for the steady state
    psolver_steady = @(State, Fluid, BC) psolver(State, Gp, Fluid, BC, rock);

    % Steady state solver
    steadystate = @(state, dt, bc) simulateToSteadyState(state, Gp, rock, fluid, dt, 'bc', bc, 'psolver', psolver_steady);
    upscale = @(Rock, BCP) upscalePermeabilityFixed(G, dp_scale, psolver, fluid_pure, Rock, L);

    % Dummy variable
    bcp = [];
end

% Do a single phase upscaling giving K which will be used for finding k_r.
perm = upscale(rock, bcp);

% Create storage for permeability and relative permeability
[kr krK] = deal(cell(2,1));
[kr{:}    krK{:}] = deal(nan(N_sat, N_dp, G.griddim^(periodic + 1)));

sat_mat = zeros(N_sat, N_dp);

lambda =@(s) bsxfun(@rdivide,fluid.relperm(s),fluid.properties());
mu = fluid.properties();

% state = initResSol(Gp, 100*barsa, 0.1);
for i = 1:N_dp
    state    = initResSol(Gp, 100*barsa, 0.0);
    dp_scale = opt.dp_vec(i);
    if periodic
        bcp.value(bcp.tags == dir) = opt.dp_vec(i);
        % Dummy variable
        bc = [];
    end
    for j = 1:N_sat
        sat = sat_vec(j);

        if ~periodic
            % If the grid is not periodic, we need to alter the bc for each
            % saturation so that the correct saturation is injected during
            % each run.
            f1 = opt.bc_faces{1};
            f2 = opt.bc_faces{2};

            fl_sat=lambda(sat*ones(G.cells.num,1));
            fl_sat=fl_sat(:,1)./sum(fl_sat,2);

            bcsat1 = fl_sat(sum(G.faces.neighbors(f1, :),2));
            bcsat2 = fl_sat(sum(G.faces.neighbors(f2, :),2));

            bc = addBC([], f1, 'pressure', dp_scale, 'sat', bcsat1);
            bc = addBC(bc, f2, 'pressure', 0       , 'sat', bcsat2);
        end


        % Calculate effective saturation for the whole domain based on fine s
        effective_sat = sum(state.s.*pv)/sum(pv(:));
        if effective_sat > 0
            tmp = (state.s - effective_sat) + repmat(sat, size(state.s));
            tmp(tmp<0) = 0;
            tmp(tmp>1) = 1;
            state.s    = tmp;
            clear tmp;
        else
            state.s = repmat(sat, size(state.s));
        end
        % Use the single phase upscaling to find timestep
        V_i = perm(dir,dir)*dp_scale/(L(dir)*min(mu));
        dt_min = 0.5*L(dir)/V_i;
        % Find steady state saturation
        [state, report] = steadystate(state, dt_min, bc);
        kr_current = fluid.relperm(state.s);

        for ph = 1:N_ph
            tmp = kr_current(:, ph);
            tmp(tmp==0)=max(tmp)*sqrt(eps);
            tmprock.perm = bsxfun(@times, rock.perm, tmp);
            tmptrans = computeTrans(G, tmprock);
            clear tmp;
            if ~all(tmptrans(Gp.cells.faces(:,end))==0)
                tmp = upscale(tmprock, bcp);
            else
                tmp = zeros(size(perm));
            end
            % We now have k_r * K. Right invert by the single phase upscaling
            % to find a value for k_r.
            krK{ph}(j,i,:) = tmp(:);
            if periodic
                % Periodic, full tensor, matrix inversion
                kr_tmp = tmp/perm;
            else
                % Diagonal, division for inversion
                kr_tmp = tmp./perm;
            end
            kr{ph}( j,i,:) = kr_tmp(:);
        end
        % There might not exist a stationary state for the desired
        % saturation value. Use the resulting value instead, weighted in
        % the usual manner by pore volume.
        sat_mat(j,i)=sum(pv.*state.s(:,1))/sum(pv);
    end
end
if N_dp == 1
    for ph = 1:N_ph
        % We remove singleton dimension if we are not really supplied a dp
        % vector. This is done to get nicer output and retain formatting
        % from an earlier implementation.
        kr{ph}  = squeeze(kr{ph});
        krK{ph} = squeeze(krK{ph});
    end
end
end

function Sp = computeMimeticIPGp(G, Gp, rock, periodic)
    Sp = computeMimeticIP(G,rock,'InnerProduct','ip_simple');
    if periodic
        hfmap = Gp.cells.faces(:,end);
        Sp.BI = Sp.BI(hfmap,hfmap);
    end
end

function Tp=computeTransGp(G,Gp,rock)
Trans=computeTrans(G,rock);
hfmap=Gp.cells.faces(:,end);
Tp=Trans(hfmap);
end


