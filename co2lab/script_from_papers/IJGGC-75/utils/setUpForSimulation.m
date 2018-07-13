function [initState, schedule, fluid, seainfo] = setUpForSimulation(Gt, rock2D, varargin)


    opt.shorter_sim_run = false;
    opt = merge_options(opt, varargin{:});
    
    fmName = 'Stofm';
    seainfo = getSeaInfo(fmName, 760 * kilo*gram/meter^3);


    % Set up fluid.
    T_ref = seainfo.seafloor_temp+273.14 + seainfo.temp_gradient * (Gt.cells.z - seainfo.seafloor_depth) / 1000; % Kelvin
    gravity on;
    p_init = seainfo.water_density * norm(gravity()) * Gt.cells.z + 1*atm;
    P_ref = mean(p_init);
    fluid = makeVEFluid(Gt, rock2D, 'sharp interface'                  , ...
                       'fixedT'       , T_ref                          , ...
                       'residual'     , [seainfo.res_sat_wat, seainfo.res_sat_co2], ...
                       'wat_rho_ref'  , seainfo.water_density          , ...
                       'co2_rho_ref'  , seainfo.rhoCref                , ... 
                       'wat_rho_pvt'  , [4.3e-5/barsa, P_ref]          , ...
                       'co2_rho_pvt'  , [[0.1, 400]*mega*Pascal, [4 250] + 274], ...
                       'wat_mu_ref'   , seainfo.water_mu               , ...
                       'co2_mu_ref'   , 0                              , ... % 0 will trigger variable viscosity
                       'co2_mu_pvt'   , [[0.1, 400]*mega*Pascal, [4 250] + 274], ...
                       'pvMult_fac'   , 1e-5/barsa                     , ...
                       'pvMult_p_ref' , P_ref                          , ...
                       'dissolution'  , false                          , ...
                       'surf_topo'    , 'smooth'                       , ...
                       'top_trap'     , []);


    % Set up initial state.
    initState.pressure = p_init;
    initState.s = repmat([1 0], Gt.cells.num, 1);
    initState.sGmax = initState.s(:,2);


    % Set up schedule, with boundary conditions and wells.
    if opt.shorter_sim_run
        itime  = 50 * year;
        isteps = 50;
        mtime1  = 10 * year;
        msteps1 = 10;
        mtime2  = 290 * year;
        msteps2 = 29;
    else
        itime  = 50 * year;
        isteps = 50;
        mtime1  = 10 * year;
        msteps1 = 10;
        mtime2  = 2900 * year;
        msteps2 = 290;
    end
    dTi         = itime / isteps;
    dTm1        = mtime1 / msteps1;
    dTm2        = mtime2 / msteps2;
    istepvec    = ones(isteps, 1) * dTi;
    mstepvec1   = ones(msteps1, 1) * dTm1;
    mstepvec2   = ones(msteps2, 1) * dTm2;
    schedule.step.val       = [istepvec; mstepvec1; mstepvec2];
    schedule.step.control   = [ones(isteps, 1); ones(msteps1+msteps2, 1) * 2];


    % Boundaries:
    % We treat all boundaries as open, however we know some boundaries are
    % fault lines that will get adjusted transmissibilities before
    % simulation
    openbfaces  = find(any(Gt.faces.neighbors==0,2));
    openbdryVal = Gt.faces.z(openbfaces) * fluid.rhoWS * norm(gravity) + 1*atm;
    bdrys = addBC( [], openbfaces, 'pressure', openbdryVal, 'sat', [1 0] );
    schedule.control(1).bc = bdrys;
    schedule.control(2).bc = bdrys;


end