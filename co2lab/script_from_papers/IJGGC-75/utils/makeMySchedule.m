function schedule = makeMySchedule(Gt, rock, fluid)
    
    surface_pressure = 0;
    
    % Fixed "injection" point, plume cells, saturations:
    pt = [4.636e5 6.4769e6]; % Utsira-South
    assert(all(size(pt)~=0));
    winx = findEnclosingCell( Gt, pt );
    

    % Setting up times for schedule:
    itime  = 30 * year;
    isteps = 30;
    mtime1  = 170 * year;
    msteps1 = 17;
    mtime2  = 2800 * year;
    msteps2 = 14;
    dTi         = itime / isteps;
    dTm1        = mtime1 / msteps1;
    dTm2        = mtime2 / msteps2;
    istepvec    = ones(isteps, 1) * dTi;
    mstepvec1   = ones(msteps1, 1) * dTm1;
    mstepvec2   = ones(msteps2, 1) * dTm2;
    schedule.step.val       = [istepvec; mstepvec1; mstepvec2];
    schedule.step.control   = [ones(isteps, 1); ones(msteps1+msteps2, 1) * 2];

    
    % Add injecting well, and shut-well:
    W = addWell([], Gt.parent, rock, winx, ...
                'name',     ['Winj' num2str(winx)],  ...
                'Type',     'rate', ...
                'Val',      1, ... % (m3/s)
                'comp_i',   [0,1], ...
                'Radius',   0.3);
    W = convertwellsVE(W, Gt.parent, Gt, rock, 'ip_tpf');
    W_shut = W;
    W_shut.val = sqrt(eps);
    schedule.control(1).W = W;
    schedule.control(2).W = W_shut;
    
    
    % Setting up boundary conditions:
    bfaces  = find(any(Gt.faces.neighbors==0,2));
    bdryVal = Gt.faces.z(bfaces) * fluid.rhoWS * norm(gravity) + surface_pressure;
    schedule.control(1).bc = addBC( [], bfaces, 'pressure', bdryVal, 'sat', [1 0] );
    schedule.control(2).bc = addBC( [], bfaces, 'pressure', bdryVal, 'sat', [1 0] );
    
end