function [output] = getTrappingCapacities(varargin)
% Computes the breakdown of structural, residual, and dissolution trapping
% 
% This is assuming CO2 is injected and fills the formation completely,
% i.e., all spill paths are utilized.
% 
% The contents of this section closely follows the implementation found in
% exploreCapacity(), which is a GUI. Here, the ranges of the input
% parameters (such as temperature gradient, etc...) are used to compute the
% range of trapping capacity (such as total and structural). Note: some
% fluid parameters are fixed to those used in exploreCapacity.

% INPUTS:

% OUTPUTS:

% see also:
%   exploreParameterRanges, exploreCapacity


Gt              = varargin{1};
rock            = varargin{2};
ta              = varargin{3};

rho_co2         = varargin{4};
rho_brine       = varargin{5};
seafloor_temp   = varargin{6};
seafloor_depth  = varargin{7};


tgrad           = varargin{8};  % degrees per kilometer
press_deviation = varargin{9};  % pressure deviation from hydrostatic (percent)
res_co2_sat     = varargin{10};
res_brine_sat   = varargin{11};
dis_max         = varargin{12}; % from CO2store



    % function handle to compute CO2 density, as function of press and temp:
    % (Get the co2 property function handles from CO2props function. Then
    % assign the density function handle to rho_co2_fun.)
    %co2 = CO2props();
    co2 = CO2props('sharp_phase_boundary', true); % 'rhofile', 'rho_demo');
    rho_co2_fun = @co2.rho;

    
    % computing pressure and temperature fields (used for estimating densities)
    gravity on;
    p_hydrostatic = rho_brine .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure 
    p = p_hydrostatic.* (1 + press_deviation/100);
    t = 273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) .* tgrad ./ 1000; % Kelvin

    
    % computing structural trap heights (H1) for each cell
    trapcells     = find(ta.traps);
    H1            = zeros(Gt.cells.num, 1);
    H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
    H1=min(H1,Gt.cells.H);
    assert(all(H1<=Gt.cells.H));


    % computing the height of the part of the column not in a structural trap
    H2 = Gt.cells.H - H1;

    
    % Computing total trapping volume in structural traps (dissolved and
    % structurally trapped
    strap_pvol_tot       = Gt.cells.volumes .* H1 .* rock.poro;
    strap_pvol_co2_plume = strap_pvol_tot .* (1 - res_brine_sat);
    strap_pvol_co2_diss  = strap_pvol_tot .* res_brine_sat .* dis_max;

    strap_mass_co2 = strap_pvol_co2_plume .* rho_co2_fun(p,t) + ...
                    strap_pvol_co2_diss  .* rho_co2;
                

    % Computing total trapping volume below structural traps (dissolved and
    % residually trapped
    btrap_pvol_tot          = Gt.cells.volumes .* H2 .* rock.poro;
    btrap_pvol_co2_residual = btrap_pvol_tot .* res_co2_sat;
    btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-res_co2_sat) .* dis_max;

    btrap_mass_co2_res = btrap_pvol_co2_residual .* rho_co2_fun(p,t);
    btrap_mass_co2_dis = btrap_pvol_co2_dissol   .* rho_co2;


    % Computing total trapping capacity per cell
    tot_trap_capa = strap_mass_co2 + btrap_mass_co2_res + btrap_mass_co2_dis;

    
    % Reporting trapping capacities:
    fprintf('Total trapping capacity: %f Gtons\n', sum(tot_trap_capa) / giga / 1e3);
    fprintf('Breakdown:\n');
    fprintf('Structural: %f Gtons\n', sum(strap_mass_co2)    / giga / 1e3);
    fprintf('Residual: %f Gtons\n', sum(btrap_mass_co2_res)  / giga / 1e3);
    fprintf('Dissolved: %f Gtons\n', sum(btrap_mass_co2_dis) / giga / 1e3);

    
    % arrays of data for output:
    output.tot_trap_capa_sum      = sum(tot_trap_capa) / giga / 1e3;
    output.strap_mass_co2_sum     = sum(strap_mass_co2) / giga / 1e3;
    output.strap_mass_co2         = strap_mass_co2;
    output.btrap_mass_co2_res_sum = sum(btrap_mass_co2_res)  / giga / 1e3;
    output.btrap_mass_co2_dis_sum = sum(btrap_mass_co2_dis) / giga / 1e3;

    
    
end




