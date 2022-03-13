function [model, schedule, initState, seainfo, other] = setupProdExModel( varargin )
    % Set up the model and schedule, using the Bjarmeland grid, open
    % boundaries, 4 CO2 injection wells, and 2 brine production wells.

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
    
    opt.surface_pressure = 1 * atm;
    opt.btype           = 'pressure';   % 'pressure' for hydrostatic or 'flux' for no-flow
    opt.itime           = 50 * year;    % injection period (years)
    opt.isteps          = 50;
    opt.mtime           = 1000 * year;  % migration period (years)
    opt.msteps          = 100;

    opt = merge_options(opt, varargin{:});

    %% Get grid and rock
    [Gt, rock2D]    = getFormationTopGrid( 'Bjarmelandfm', 4 );

    
    %% Set-up fluid structure
    seainfo    = getSeaInfo('Bjarmelandfm', 760);
    gravity on;
    caprock_pressure = (Gt.cells.z * seainfo.water_density * norm(gravity)) ...
        .* (1 + seainfo.press_deviation/100) + opt.surface_pressure;
    ta = trapAnalysis(Gt, true);
    co2 = CO2props();
    ref_temp = seainfo.seafloor_temp + 273.15;
    T_ref = ref_temp + seainfo.temp_gradient * (Gt.cells.z - seainfo.seafloor_depth) / 1000;   
    fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                               , ...
                       'fixedT'       , T_ref                                       , ...
                       'residual'     , [seainfo.res_sat_wat , seainfo.res_sat_co2] , ...
                       'wat_rho_ref'  , seainfo.water_density                       , ...
                       'co2_rho_ref'  , seainfo.rhoCref                             , ...
                       'wat_rho_pvt'  , [4.3e-5/barsa   , mean(caprock_pressure)]   , ...
                       'co2_rho_pvt'  , [[0.1, 400] * mega * Pascal   , [4 250] + 274]  , ...
                       'wat_mu_ref'   , seainfo.water_mu                                , ...
                       'co2_mu_pvt'   , [[0.1, 400] * mega * Pascal   , [4 250] + 274]  , ...
                       'pvMult_fac'   , 1e-5/barsa                                  , ...
                       'pvMult_p_ref' , mean(caprock_pressure)                      , ...
                       'dissolution'  , false                                       , ...
                       'dis_rate'     , []                                          , ...
                       'dis_max'      , []                                          , ...
                       'surf_topo'    , 'smooth'                                    , ...
                       'top_trap'     , []);

                   
    %% Set-up initial state
    initState.pressure = fluid.rhoWS * norm(gravity()) * Gt.cells.z + opt.surface_pressure;
    initState.s = repmat([1 0], Gt.cells.num, 1);
    initState.sGmax = initState.s(:,2);

    
    %% Set-up wells
    [P_over] = computeOverburdenPressure(Gt, rock2D, seainfo.seafloor_depth, fluid.rhoWS);

    % Wells:

    %   ---get cell indexes from physical coordinate:
    inj_wcoord  = [ 1.0170e6 8.0883e6; 1.0470e6 8.1703e6; 1.0750e6 8.1743e6; 1.0690e6 8.1323e6 ]; % physical coordinates
    prd_wcoord  = [ 1.0570e6 8.1963e6; 1.0430e6 8.1043e6 ];
    inj_wc = findEnclosingCell(Gt, inj_wcoord);
    prd_wc = findEnclosingCell(Gt, prd_wcoord);

    %   ---explicitly set injection rate and production bhp:
    inj_rate    = [0.0515 0.0088 0.0661 0.0887]; % m3/s
    prd_bhp     = [initState.pressure(prd_wc(1))/2, initState.pressure(prd_wc(2))/2]; % Pascals

    assert( size(inj_wcoord,1) == numel(inj_rate) )
    assert( size(prd_wcoord,1) == numel(prd_bhp) )

    %   ---create well structure: (could use addWellVE)
    W = [];
    for i=1:4
       W = addWell(W, Gt.parent, rock2D, inj_wc(i), ...
           'name',      sprintf('Winj%i', inj_wc(i)), ...
           'sign',      1, ...
           'Type',      'rate', ...
           'Val',       inj_rate(i), ...
           'comp_i',    [0 1], ...
           'Radius',    0.3 );
    end
    for i=1:2
       W = addWell(W, Gt.parent, rock2D, prd_wc(i), ...
           'name',      sprintf('Wprd%i', prd_wc(i)), ...
           'sign',      -1, ...
           'Type',      'bhp', ...
           'Val',       prd_bhp(i), ...
           'comp_i',    [1 0], ...
           'Radius',    0.3 );
    end
    W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');

    % Make adjustments to well controls for migration period
    W_shut = W;
    for i=1:4
        W_shut(i).val = sqrt(eps);
    end
    for i=5:6
        W_shut(i).type = 'rate';
        W_shut(i).val  = -sqrt(eps);
    end

    
    %% Construct schedule
    schedule.control(1).W = W;
    schedule.control(2).W = W_shut;
    dTi         = opt.itime / opt.isteps;
    dTm         = opt.mtime / opt.msteps;
    istepvec    = ones(opt.isteps, 1) * dTi;
    mstepvec    = ones(opt.msteps, 1) * dTm;
    schedule.step.val       = [istepvec; mstepvec];
    schedule.step.control   = [ones(opt.isteps, 1); ones(opt.msteps, 1) * 2];

    
    %% Add boundary conditions (either constant pressure, or no-flow)
    btype  = opt.btype;
    bfaces = find(any(Gt.faces.neighbors==0,2)); 
    if strcmpi(btype,'pressure')
        bdryVal     = Gt.faces.z(bfaces) * fluid.rhoWS * norm(gravity) + opt.surface_pressure;
    elseif strcmpi(btype,'flux')
        bdryVal     = zeros(numel(bfaces),1);
    end
    bc = addBC( [], bfaces, btype, bdryVal, 'sat', [1 0] );
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end
    
    
    %% Make model
    model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

    
    %% Possibly required variables (especially for post-processing):
    other.ta                = ta;
    other.surface_pressure  = opt.surface_pressure;
    other.btype             = opt.btype;
    other.rock = rock2D;
    other.residual = [fluid.res_water, fluid.res_gas];
    other.ref_temp = seainfo.seafloor_temp + 273.15;
    other.ref_depth = seainfo.seafloor_depth;
    other.temp_grad = seainfo.temp_gradient;
    other.rhoW = seainfo.water_density;

end
