function [ model, schedule, initState, seainfo, other ] = makeBjarmelandModel( varargin )
% Simple Bjarmeland model for testing pressure-limited injection
% for CO2 storage.
%
% Note: The Bjarmeland model is part of the Norwegian Petroleum
% Directorate's Compiled CO2 Storage Atlas of the Norwegian Continental
% Shelf, found here:
%       http://www.npd.no/en/Publications/Reports/Compiled-CO2-atlas/

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

    opt.coarsening      = 10;
    opt.surface_pressure = 1 * atm;
    opt.btype           = 'pressure';   % 'pressure' for hydrostatic or 'flux' for no-flow
    opt.itime           = 30;           % injection period (years)
    opt.isteps          = 30;
    opt.mtime           = 500;         % migration period (years)
    opt.msteps          = 50;
    opt.wcoord          = [1.0250e6 8.0543e6; ...
                           1.0590e6 8.1543e6];      % well coordinates [x1 y1; x2 y2; ...]
    opt.qt              = [3.7028e12; 2.6198e12];   % total injected mass per well (kg of CO2)
                                                    % can be taken from the structural trapping 
                                                    % capacity upstream from well location, or 
                                                    % a scaled value
    opt.wvals = []; % if empty, uses opt.qt, otherwise uses opt.wvals
    
    opt = merge_options(opt, varargin{:});
    assert( size(opt.wcoord,1) == numel(opt.qt) )
    gravity on;
    
    %% Construct model
    % Grid and rock:
    [Gt, rock2D]        = getFormationTopGrid('Bjarmelandfm', opt.coarsening );
    dh = [];
    
    % Set up fluid:
    rhoCref             = 760 * kilogram / meter ^3;
    seainfo             = getSeaInfo('Bjarmelandfm', rhoCref);
    caprock_pressure    = (Gt.cells.z * seainfo.water_density * norm(gravity)) ...
                            .* (1 + seainfo.press_deviation/100) + opt.surface_pressure;  
    T_ref               = seainfo.seafloor_temp + 273.15 + seainfo.temp_gradient ...
                            * (Gt.cells.z - seainfo.seafloor_depth) / 1000; % Kelvin
    fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                           , ...
                        'fixedT'       , T_ref                                  , ...
                        'residual'     , [seainfo.res_sat_co2, seainfo.res_sat_wat] , ...
                        'wat_rho_ref'  , seainfo.water_density                  , ...
                        'co2_rho_ref'  , seainfo.rhoCref                        , ...
                        'wat_rho_pvt'  , [4.3e-5/barsa, mean(caprock_pressure)] , ...
                        'co2_rho_pvt'  , [[0.1, 400]*mega*Pascal, [4 250]+274]  , ...
                        'wat_mu_ref'   , seainfo.water_mu                       , ...
                        'co2_mu_pvt'   , [[0.1, 400]*mega*Pascal, [4 250]+274]  , ...
                        'pvMult_fac'   , 1e-5/barsa                             , ...
                        'pvMult_p_ref' , mean(caprock_pressure)                 , ... 
                        'dissolution'  , false                                  , ...
                        'dis_rate'     , 0                                      , ... % 0 means instantaneous
                        'dis_max'      , 53/seainfo.rhoCref                     , ... % i.e., 53/760 = 0.07; % 1 kg water holds 0.07 kg of CO2
                        'surf_topo'    , 'smooth'                               , ...
                        'top_trap'     , dh);
             
    model       = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
    
    % possibly required variables (especially for post-processing):
    other.ta                = trapAnalysis(Gt, true);
    other.surface_pressure  = opt.surface_pressure;
    other.btype             = opt.btype;
    other.rock = rock2D;
    other.dh = dh;
    other.residual = [fluid.res_water, fluid.res_gas];
    other.ref_temp = seainfo.seafloor_temp + 273.15;
    other.ref_depth = seainfo.seafloor_depth;
    other.temp_grad = seainfo.temp_gradient;
    other.rhoW = seainfo.water_density;
                   
    %% Set initial conditions:
    %   (hydrostatic pressure, water saturated)
    initState.pressure  = fluid.rhoWS * norm(gravity()) * Gt.cells.z + opt.surface_pressure;
    initState.s         = repmat([1 0], Gt.cells.num, 1);
    initState.sGmax     = initState.s(:,2);
    
    %% Set up schedule:
    %   including injection wells & rates, migration rates, time steps, boundary conditions
    wqtots  = opt.qt./fluid.rhoGS;                  % this could be scaled down. % total volume (m3) to inject per well in terms of reference density
    if isempty(opt.wvals)
        wvals   = wqtots./convertFrom(opt.itime,year);  % injection rates (m3/s) with respect to CO2 reference density
    else
        wvals = opt.wvals;
    end
    wcells = findEnclosingCell(Gt, opt.wcoord);
    assert( isfield(Gt,'parent') )
    W = [];
    for i=1:numel(wvals)
        W = addWell(W, Gt.parent, rock2D, wcells(i), ...
                    'name',     sprintf('Winj%i', wcells(i)),  ...
                    'Type',     'rate', ...
                    'Val',      wvals(i), ...
                    'comp_i',   [0 1], ...
                    'Radius',   0.3 );
    end
    W = convertwellsVE(W, Gt.parent, Gt, rock2D, 'ip_tpf');
    W_shut = W;
    for i = 1:numel(W_shut)
        W_shut(i).type  = 'rate';
        W_shut(i).val   = sqrt(eps); % the minimum value
    end
    schedule.control(1).W = W;
    schedule.control(2).W = W_shut;
    dTi         = convertFrom(opt.itime,year) / opt.isteps;
    dTm         = convertFrom(opt.mtime,year) / opt.msteps;
    istepvec    = ones(opt.isteps, 1) * dTi;
    mstepvec    = ones(opt.msteps, 1) * dTm;
    schedule.step.val       = [istepvec; mstepvec];
    schedule.step.control   = [ones(opt.isteps, 1); ones(opt.msteps, 1) * 2];
    %%{
    % Optional plotting to check placement and rates:
    figure(100); clf; set(gcf,'Position',[1 1 953 615])
    subplot(2,2,[2 4])
    ta = other.ta;
    mapPlot(gcf, Gt, 'traps', ta.traps, ...
        'trapcolor', [0.5 0.5 0.5], 'trapalpha', 0.7, ...
        'rivers', ta.cell_lines, 'rivercolor', [0 0 1], ...
        'maplines', 40, 'wellcells', wcells, 'well_numbering', true);
    colorizeCatchmentRegions(Gt, ta);
    plotGrid(Gt, 'FaceColor','none','EdgeAlpha',0.1)
    axis equal tight off
    rates = [schedule.control(1).W.val].*fluid.rhoGS.*(1*year)./10^9; % Mt/yr
    subplot(2,2,[1 3]); bar(rates);
    xlim([0 numel(rates)+1]);
    xlabel('well number');
    ylabel({'initial rate [Mt/yr]';['for ',num2str(opt.itime,year),' yrs inj.']});
    %}
    
    % Boundary conditions:
    %   (either constant pressure, or no-flow)
    bfaces = find(any(Gt.faces.neighbors==0,2));
    if strcmpi(opt.btype,'pressure')
        bdryVal     = Gt.faces.z(bfaces) * fluid.rhoWS * norm(gravity) + opt.surface_pressure;
    elseif strcmpi(opt.btype,'flux')
        bdryVal     = zeros(numel(bfaces),1);
    end
    bc = addBC( [], bfaces, opt.btype, bdryVal, 'sat', [1 0] );
    for i = 1:numel(schedule.control)
        schedule.control(i).bc = bc;
    end
end
