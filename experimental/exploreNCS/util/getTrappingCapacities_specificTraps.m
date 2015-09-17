function [output] = getTrappingCapacities_specificTraps(varargin)
% Computes the trapping details for a specific structural trap, which is
% identified by the trapping ID passed into function. The details computed
% include:
%   rock volume of structural trap
%   net volume (based on NTG)
%   pore volume (based on porosity and net volume)
%   average ntg, perm, poro and trap depth (min and max)
%   storage capacity of trap, assuming 100% storage efficiency, in mass and
%   volume

% INPUTS:

% OUTPUTS:

% see also:
%   exploreParameterRanges, exploreCapacity, getTrappingCapacities


Gt              = varargin{1};
rock            = varargin{2};
ta              = varargin{3};

trapID          = varargin{4};
trapName        = varargin{5}; % needs to be a cell

rho_co2         = varargin{6};
rho_brine       = varargin{7};
seafloor_temp   = varargin{8};
seafloor_depth  = varargin{9};


tgrad           = varargin{10};  % degrees per kilometer
press_deviation = varargin{11};  % pressure deviation from hydrostatic (percent)

res_co2_sat     = varargin{12};
res_brine_sat   = varargin{13};
dis_max         = varargin{14}; % from CO2store



    % function handle to compute CO2 density, as function of press and temp:
    % (Get the co2 property function handles from CO2props function. Then
    % assign the density function handle to rho_co2_fun.)
    %co2 = CO2props();
    co2 = CO2props('sharp_phase_boundary', true, 'rhofile', 'rho_demo');
    rho_co2_fun = @co2.rho;

    
    % computing pressure and temperature fields (used for estimating densities)
    gravity on;
    p_hydrostatic = rho_brine .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure 
    p = p_hydrostatic.* (1 + press_deviation/100);
    t = 273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) .* tgrad ./ 1000; % Kelvin

    
    % computing structural trap heights (H1) for each cell in the specified
    % structural trap (by the trapID)
    trapcells     = find(ta.traps == trapID);
    H1            = zeros(Gt.cells.num, 1);
    H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
    H1=min(H1,Gt.cells.H);
    assert(all(H1<=Gt.cells.H));


%     % computing the height of the part of the column not in the specified
%     % structural trap
%     H2 = Gt.cells.H - H1;
% 
%     
%     % Computing total trapping volume in the specified structural trap
%     % (dissolved and structurally trapped)
%     % Note: rock.poro is an array of values
%     strap_pvol_tot       = Gt.cells.volumes .* H1 .* rock.poro; % cubic meters, each cell in trap
%     strap_pvol_co2_plume = strap_pvol_tot .* (1 - res_brine_sat);
%     strap_pvol_co2_diss  = strap_pvol_tot .* res_brine_sat .* dis_max;
% 
%     strap_mass_co2 = strap_pvol_co2_plume .* rho_co2_fun(p,t) + ...
%                     strap_pvol_co2_diss  .* rho_co2;
                
                
    % Other info about The Structural Trap that might be useful for
    % comparions later:
    
    % number of cells:
    output.numtrapcells = numel(trapcells);
    
    % rock vol:
    output.rockVol_m3 = sum(Gt.cells.volumes .* H1);
    
    % net vol:
    output.netVol_m3 = sum(Gt.cells.volumes .* H1 .* rock.ntg);
    
    % pore vol: (included NTG)
    output.poreVol_m3 = sum(Gt.cells.volumes .* H1 .* rock.ntg .* rock.poro); % cubic meters, total in structural trap
    
    % av ntg, poro, perm, and their ranges:
    output.ntg_min = min(rock.ntg(trapcells));
    output.ntg_max = max(rock.ntg(trapcells));
    output.ntg_av  = mean(rock.ntg(trapcells));
    
    output.perm_min = min(rock.perm(trapcells));
    output.perm_max = max(rock.perm(trapcells));
    output.perm_av  = mean(rock.perm(trapcells));
    
    output.poro_min = min(rock.poro(trapcells));
    output.poro_max = max(rock.poro(trapcells));
    output.poro_av  = mean(rock.poro(trapcells));
    
    % min/max depth of structural trap (not formation), and average:
    output.depth_min = min(Gt.cells.z(trapcells));
    output.depth_max = max(Gt.cells.z(trapcells));
    output.depth_av  = mean(Gt.cells.z(trapcells));
    

    % storage capacity in structure, in volume (m3), and in mass (kg, Mt)
    % assuming 100% of pore volume is capable of storing CO2 (no storage
    % efficiency factor used here):
    output.storageCap_m3 = output.poreVol_m3;
    output.storageCap_kg = sum( Gt.cells.volumes(trapcells) .* H1(trapcells) .* rock.ntg(trapcells) .* rock.poro(trapcells) .* rho_co2_fun(p(trapcells),t(trapcells)) ); % m3 * kg/m3 = kg
    output.storageCap_Mt = output.storageCap_kg / 1e9; % 10^9 kg = 1 Mt
    
    
    % corresponding storage efficiency in structural trap:
    output.Seff = output.storageCap_m3 / output.poreVol_m3;
    
    
    
    output.T = table(  [output.numtrapcells;        ...
                output.rockVol_m3;          ...
                output.netVol_m3;           ...
                output.poreVol_m3;          ...
                output.ntg_min;             ...
                output.ntg_max;             ...
                output.ntg_av;              ...
                output.perm_min;            ...
                output.perm_max;            ...
                output.perm_av;             ...
                output.poro_min;            ...
                output.poro_max;            ...
                output.poro_av;            ...
                output.depth_min;           ...
                output.depth_max;           ...
                output.depth_av;            ...
                output.storageCap_Mt],       ...
        'RowNames', {'NumTrapCells' 'RockVolm3' 'NetVolm3' 'PoreVolm3' ...
        'ntgmin' 'ntgmax' 'ntgav' 'permmin' 'permmax' 'permav' 'poromin' 'poromax' 'poroav' ...
        'depthmin' 'depthmax' 'depthav' 'StorageCapacityMt'}, ...
        'VariableNames', trapName  );
    
    output.trapName = trapName;
    
    
%     T = table(  output.numtrapcells,        ...
%                 output.rockVol_m3,          ...
%                 output.netVol_m3,           ...
%                 output.poreVol_m3,          ...
%                 output.ntg_min,             ...
%                 output.ntg_max,             ...
%                 output.ntg_av,              ...
%                 output.perm_min,            ...
%                 output.perm_max,            ...
%                 output.perm_av,             ...
%                 output.poro_min,            ...
%                 output.poro_max,            ...
%                 output.poro_av,             ...
%                 output.depth_min,           ...
%                 output.depth_max,           ...
%                 output.depth_av,            ...
%                 output.storageCap_Mt,       ...
%         'RowNames', trapName, ...
%         'VariableNames', {'NumTrapCells' 'RockVolm3' 'NetVolm3' 'PoreVolm3' ...
%         'ntgmin' 'ntgmax' 'ntgav' 'permmin' 'permmax' 'permav' 'poromin' 'poromax' 'poroav' ...
%         'depthmin' 'depthmax' 'depthav' 'StorageCapacityMt'} )
%     
    
    

%     % Computing total trapping volume below the specified structural trap
%     % (dissolved and residually trapped
%     btrap_pvol_tot          = Gt.cells.volumes .* H2 .* rock.poro;
%     btrap_pvol_co2_residual = btrap_pvol_tot .* res_co2_sat;
%     btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-res_co2_sat) .* dis_max;
% 
%     btrap_mass_co2_res = btrap_pvol_co2_residual .* rho_co2_fun(p,t);
%     btrap_mass_co2_dis = btrap_pvol_co2_dissol   .* rho_co2;
% 
% 
%     % Computing total trapping capacity per cell
%     tot_trap_capa = strap_mass_co2 + btrap_mass_co2_res + btrap_mass_co2_dis;
% 
%     
%     % Reporting trapping capacities:
%     fprintf('Total trapping capacity: %f Gtons\n', sum(tot_trap_capa) / giga / 1e3);
%     fprintf('Breakdown:\n');
%     fprintf('Structural: %f Gtons\n', sum(strap_mass_co2)    / giga / 1e3);
%     fprintf('Residual: %f Gtons\n', sum(btrap_mass_co2_res)  / giga / 1e3);
%     fprintf('Dissolved: %f Gtons\n', sum(btrap_mass_co2_dis) / giga / 1e3);
% 
%     
%     % arrays of data for output:
%     output.tot_trap_capa_sum      = sum(tot_trap_capa) / giga / 1e3;
%     output.strap_mass_co2_sum     = sum(strap_mass_co2) / giga / 1e3;
%     output.strap_mass_co2         = strap_mass_co2;
%     output.btrap_mass_co2_res_sum = sum(btrap_mass_co2_res)  / giga / 1e3;
%     output.btrap_mass_co2_dis_sum = sum(btrap_mass_co2_dis) / giga / 1e3;

    
    
end




