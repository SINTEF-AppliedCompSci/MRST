function [sat_vec, kr, krK, perm, varargout] = ...
          upscaleRelpermLimit(G, rock, fluid, varargin)
% Upscale relperm based on viscous/capillary limit tabulated by saturation
%
% SYNOPSIS:
%   calculateRelperm(G, rock, fluid)
%   calculateRelperm(G, rock, fluid, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure as described by grid_structure, or
%             alternatively a periodic grid as created by
%             makePeriodicGridMulti3d. If using a periodic grid the
%             keyword 'type' must be set to 'periodic'
%
%   rock    - Valid rock data structure.
%
%   fluid   - Fluid object as defined by initSWOFFluidJfunc.
%
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%           - type -- The type of grid provided. Can be either 'fixed' for
%                     a regular grid or 'periodic' for a periodic grid.
%                     Defaults to 'fixed.
%
%           - n    -- The number of saturation values for which the relperm
%                     is to be calculated. Defaults to 100, divides the
%                     interval [0,1] into N equal increments.
%
%           - limit -- Limiting value for the calculation. 'capillary'
%                      capillary upscaling, 'viscious' gives viscosity
%                      based upscaling.
%
% RETURNS:
%   sat_vec  - The saturation values used for upscaling. Useful for
%              plotting and tabulation.
%
%  kr        - Relative permeability values
%
%  perm    - Permeability upscaled using pure fluid.
%
%  krK       - Effective permeability. The upscaling is based on first
%              finding an upscaled permeability with a fluid with unitary
%              values (fluid_pure). Once this is found, relperm is found by
%              solving the expression kr * K = krK. Depending on the
%              permeability type (full tensor / diagonal) this is be done
%              by inversion of K or division.
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


opt = struct('type' ,   'fixed', ...
             'n'    ,   100, ...
             'limit',   'capillary'...
             );

opt         = merge_options(opt, varargin{:});
d           = G.griddim;
dp_scale    = 1e-3;
periodic    = strcmpi(opt.type, 'periodic');
L           = max(G.faces.centroids)-min(G.faces.centroids);
fluid_pure  = initSingleFluid('mu' ,1, 'rho', 1);
pv          = poreVolume(G,rock);

% Pairs of opposing boundary condition keywords
bcsides = {'LEFT', 'RIGHT';...
           'BACK', 'FRONT';...
           'TOP', 'BOTTOM'};

IP      = @(Grid, Rock) computeMimeticIPGp(G, Grid, Rock, periodic);
psolver = @(state, Grid, Fluid, BCP, BC, Rock) ...
   incompMimetic(state, Grid, IP(Grid, Rock), Fluid, 'bcp', BCP, 'bc', BC);

if periodic
    bcl   = cell(d,1);
    bcr   = cell(d,1);
    dp    = cell(d,1);

    for i = 1:d
        bcl{i} = pside([], G, bcsides{i, 1}, 0);
        bcr{i} = pside([], G, bcsides{i, 2}, 0);
        dp{i} = 0;
    end

    dp{1} = dp_scale;
    % Create a periodic grid which couples the edges so that bcl{i}.faces
    % are coupled with bcr{i}.faces for i \in G.griddim

    [Gp, bcp] = makePeriodicGridMulti3d(G, bcl, bcr, dp);

    % Define a pressure solver based on the earlier function which uses
    % periodic boundary conditions
    psolver = @(State, Grid, Fluid, BCP, Rock) psolver(State, Grid, Fluid, BCP, [], Rock);
    upscale = @(Rock) upscalePermeabilityPeriodic(Gp, bcp, dp_scale, psolver, fluid_pure, Rock, L);
    clear bcl bcr
else
    % Fixed grid, will not create a periodic grid.
    Gp = G;
    % Define pressure solver which uses normal BC
    psolver = @(State, Grid, Fluid, BC, Rock) psolver(State, Grid, Fluid, [], BC, Rock);
    upscale = @(Rock) upscalePermeabilityFixed(G, dp_scale, psolver, fluid_pure, Rock, L);
end

perm = upscale(rock);
if strcmpi(opt.limit, 'viscous')
    [sat_vec, kr, krK] = upscaleViscous(Gp, G, opt.n, periodic, pv, fluid, upscale, rock, perm);
else
    [sat_vec, kr, krK, pc_values] = upscaleCapillary(Gp, G, opt.n, periodic, pv, fluid, upscale, rock, perm);
    varargout{1} = pc_values;
end

end

function Sp = computeMimeticIPGp(G, Gp, rock, periodic)
    Sp = computeMimeticIP(G,rock,'InnerProduct','ip_simple');
    if periodic
        hfmap = Gp.cells.faces(:,end);
        Sp.BI = Sp.BI(hfmap,hfmap);
    end
end

function [sat_total, kr, krK, pc_values] = upscaleCapillary(G, G_parent, N, periodic, pv, fluid, upscale, rock, perm)

    [pc_upsc, pc_max, pc_min, sat_min, sat_max] = upscalePCCaplimit(G_parent, fluid, pv);
    sat_points = linspace(sat_min, sat_max, N);
    pc_values  = pc_upsc(sat_points);
    ix = ~isnan(pc_values);
    pc_values = pc_values(ix);
    sat_points = sat_points(ix);
    sat_total=nan(N,1);

    % Helper column of size G.cells.num for faking saturations
    col = ones(G.cells.num,1);

    % Create storage for permeability and relative permeability
    [kr krK] = deal(cell(2,1));
    [kr{:}    krK{:}] = deal(nan(N, G.griddim^(periodic + 1)));


    for i=1:numel(sat_points)
       sat=fluid.pcinv(pc_values(i).*col);
       if(all(sat<=1 & sat>=0))
          sat_total(i)=sum(pv.*sat)/sum(pv);
          kr_tmp=fluid.relperm(sat);

          for kk=1:2;
             % hack in case of non connected domains of fluid
             % should be fixed in pressure solver
             max_kr = max(kr_tmp(:,kk));
             tmp     = max(max_kr*1e-10, kr_tmp(:,kk));

             % multilply permeability with relperm
             tmprock.perm=bsxfun(@times,rock.perm,tmp);
             clear tmp;

             %only done to test to main
             T = computeTrans(G_parent, tmprock);
             TT=T(G_parent.cells.faces(:,end));

             % do single phase upscaling in all direction given the stationary
             % state
             if(~all(TT==0));
                Kkr_tmp = upscale(tmprock);
             else
                Kkr_tmp = zeros(G.griddim, G.griddim^(periodic));
             end
             krK{kk}(i,:)=Kkr_tmp(:);

             if periodic
                tmp = Kkr_tmp/perm;
             else
                tmp = Kkr_tmp./perm;
             end
             kr{kk}(i,:)=tmp(:);
             clear tmp;
          end
       end
    end
end

function [sat_total, kr, krK] = upscaleViscous(G, G_parent, N, periodic, pv, fluid, upscale, rock, perm)
    fl_vec    = linspace(0,1,N);
    fl_inv    = @(f) fluid.f_inv(f);
    sat_total = nan(N,1);

    % Helper column of size G.cells.num for faking saturations
    col = ones(G.cells.num,1);

    % Create storage for permeability and relative permeability
    [kr krK] = deal(cell(2,1));
    [kr{:}    krK{:}] = deal(nan(N, G.griddim^(periodic + 1)));


    mu = fluid.properties();

    for i=1:N
       % Calculate mobility based on current tabulated saturation
       tmp = fluid.relperm(fl_vec(i).*col);
       lam = bsxfun(@rdivide, tmp, mu);
       clear tmp;

       f   = sum(pv.*(lam(:,1)./sum(lam, 2)))/sum(pv);
       sat = fl_inv(f.*col);

       if all(sat<= 1 & sat >= 0)
          sat_total(i) = sum(pv.*sat)/sum(pv);
          kr_tmp       = fluid.relperm(sat);

          for kk=1:2;
             % hack in case of non connected domains of fluid
             % should be fixed in pressure solver
             max_kr = max(kr_tmp(:,kk));
             tmp     = max(max_kr*1e-10, kr_tmp(:,kk));

             % multilply permeability with relperm
             tmprock.perm=bsxfun(@times,rock.perm,tmp);

             % only done to test to main
             T  = computeTrans(G_parent,tmprock);
             TT = T(G_parent.cells.faces(:,end));
             % do single phase upscaling in all direction given the stationary
             % state
             if(~all(TT==0));
                Kkr_tmp = upscale(tmprock);
             else
                Kkr_tmp = zeros(G.griddim, G.griddim^(periodic));
             end

             krK{kk}(i,:) = Kkr_tmp(:);
             if periodic
                tmp = Kkr_tmp/perm;
             else
                tmp = Kkr_tmp./perm;
             end
             kr{kk}(i,:) = tmp(:);
             clear tmp tmprock
          end
       end
    end
end
