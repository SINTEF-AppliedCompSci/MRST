function [output] = exploreParameterRanges(varargin)
% Computes the breakdown of structural, residual, and dissolution trapping.
%
% Ex: output = exploreParameterRanges(Gt, rock, ta, ...
%                       parameter2vary, parameterRange)
%
% parameter2vary and parameterRange are optional. If not specified, no
% range is used but rather parameters are set to their default values.
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
%   parameter2vary - either 'tgrad', 'press_deviation',
%   'res_co2_sat', 'res_brine_sat', or 'dis_max'

% OUTPUTS:



Gt = varargin{1};
rock = varargin{2};
ta = varargin{3};



% Fixed Inputs: (co2_rho, rho_brine, poro, seafloor temp and depth)
rho_co2         = 760 * kilogram / meter ^3; % reference density (while co2 density is function of p,t)
rho_brine       = 1000; % kg per m3
seafloor_temp   = 7; % in Celsius
seafloor_depth  = 100 * meter; % meter % in North Sea?


% Default values:
tgrad_default           = 35.6;         % degrees per kilometer
press_deviation_default = 0;            % pressure deviation from hydrostatic (percent)
res_co2_sat_default     = 0.21;
res_brine_sat_default   = 0.11;
diss_max_default        = (53 * kilogram / meter^3) / rho_co2; % from CO2store   



if numel(varargin)==5
    parameter2vary = varargin{4};
    parameterRange = varargin{5};

    % ************ Inputs to Vary: ************
    % Set inputs: (will use either _default or _range)
    if strcmpi(parameter2vary,'tgrad')
        tgrad_input = parameterRange; numPts = numel(tgrad_input);

        % tgrad
        press_deviation_input(1:numPts) = press_deviation_default;
        res_co2_sat_input(1:numPts)     = res_co2_sat_default;
        res_brine_sat_input(1:numPts)   = res_brine_sat_default;
        dis_max_input(1:numPts)         = diss_max_default;

    elseif strcmpi(parameter2vary,'press_deviation')
        press_deviation_input = parameterRange; numPts = numel(press_deviation_input);

        tgrad_input(1:numPts)           = tgrad_default;
        % press
        res_co2_sat_input(1:numPts)     = res_co2_sat_default;
        res_brine_sat_input(1:numPts)   = res_brine_sat_default;
        dis_max_input(1:numPts)         = diss_max_default;

    elseif strcmpi(parameter2vary,'res_co2_sat')
        res_co2_sat_input  = parameterRange; numPts = numel(res_co2_sat_input);

        tgrad_input(1:numPts)           = tgrad_default;
        press_deviation_input(1:numPts) = press_deviation_default;
        % res_co2
        res_brine_sat_input(1:numPts)   = res_brine_sat_default;
        dis_max_input(1:numPts)         = diss_max_default;

    elseif strcmpi(parameter2vary,'res_brine_sat')
        res_brine_sat_input = parameterRange; numPts = numel(res_brine_sat_input);

        tgrad_input(1:numPts)           = tgrad_default;
        press_deviation_input(1:numPts) = press_deviation_default;
        res_co2_sat_input(1:numPts)     = res_co2_sat_default;
        % res_brine
        dis_max_input(1:numPts)         = diss_max_default;

    elseif strcmpi(parameter2vary,'dis_max')
        dis_max_input = parameterRange * (kilogram/meter^3)/rho_co2; numPts = numel(dis_max_input);

        tgrad_input(1:numPts)           = tgrad_default;
        press_deviation_input(1:numPts) = press_deviation_default;
        res_co2_sat_input(1:numPts)     = res_co2_sat_default;
        res_brine_sat_input(1:numPts)   = res_brine_sat_default;
        % dis_max
        
    end
    % *****************************************

else
    disp('Parameter2vary input is undefined: using default values instead.')

    tgrad_input             = tgrad_default;
    press_deviation_input   = press_deviation_default;
    res_co2_sat_input       = res_co2_sat_default;
    res_brine_sat_input     = res_brine_sat_default;
    dis_max_input           = diss_max_default;
    numPts                  = 1;

end
    
    
    



tot_trap_capa_sum       = zeros(numPts,1);
strap_mass_co2_sum      = zeros(numPts,1);
btrap_mass_co2_res_sum  = zeros(numPts,1);
btrap_mass_co2_dis_sum  = zeros(numPts,1);

for i = 1:numPts;

    tgrad           = tgrad_input(i);
    press_deviation = press_deviation_input(i);
    res_co2_sat     = res_co2_sat_input(i);
    res_brine_sat   = res_brine_sat_input(i);
    dis_max         = dis_max_input(i);

    % function handle to compute CO2 density, as function of press and temp:
    % (Get the co2 property function handles from CO2props function. Then
    % assign the density function handle to rho_co2_fun.)
    %co2 = CO2props();
    co2 = CO2props('sharp_phase_boundary', true, 'rhofile', 'rho_demo');
    rho_co2_fun = @co2.rho;

    % computing pressure and temperature fields (used for estimating densities)
    %fprintf('Computing pressure and temperature ... ');
    gravity on;
    p_hydrostatic = rho_brine .* norm(gravity) .* Gt.cells.z; % hydrostatic pressure 
    p = p_hydrostatic.* (1 + press_deviation/100);
    t = 273.15 + seafloor_temp + (Gt.cells.z - seafloor_depth) .* tgrad ./ 1000; % Kelvin
    %fprintf('done\n');

    % computing structural trap heights (H1) for each cell
    %fprintf('Computing structural trap heights ... ');
    trapcells     = find(ta.traps);
    H1            = zeros(Gt.cells.num, 1);
    H1(trapcells) = ta.trap_z(ta.traps(trapcells)) - Gt.cells.z(trapcells);
    H1=min(H1,Gt.cells.H);
    assert(all(H1<=Gt.cells.H));
    %fprintf('done\n');

    % computing the height of the part of the column not in a structural trap
    H2 = Gt.cells.H - H1;

    % Computing total trapping volume in structural traps (dissolved and
    % structurally trapped
    %fprintf('Computing total volume in structural traps ... ');
    strap_pvol_tot       = Gt.cells.volumes .* H1 .* rock.poro;
    strap_pvol_co2_plume = strap_pvol_tot .* (1 - res_brine_sat);
    strap_pvol_co2_diss  = strap_pvol_tot .* res_brine_sat .* dis_max;

    strap_mass_co2 = strap_pvol_co2_plume .* rho_co2_fun(p,t) + ...
                    strap_pvol_co2_diss  .* rho_co2;
    %fprintf('done\n');

    % Computing total trapping volume below structural traps (dissolved and
    % residually trapped
    %fprintf('Computing trapping volume below structural traps ... ');
    btrap_pvol_tot          = Gt.cells.volumes .* H2 .* rock.poro;
    btrap_pvol_co2_residual = btrap_pvol_tot .* res_co2_sat;
    btrap_pvol_co2_dissol   = btrap_pvol_tot .* (1-res_co2_sat) .* dis_max;

    btrap_mass_co2_res = btrap_pvol_co2_residual .* rho_co2_fun(p,t);
    btrap_mass_co2_dis = btrap_pvol_co2_dissol   .* rho_co2;
    %fprintf('done\n');

    % Computing total trapping capacity per cell
    tot_trap_capa = strap_mass_co2 + btrap_mass_co2_res + btrap_mass_co2_dis;

    % Reporting trapping capacities:
    fprintf('Total trapping capacity: %f Gtons\n', sum(tot_trap_capa) / giga / 1e3);
    fprintf('Breakdown:\n');
    fprintf('Structural: %f Gtons\n', sum(strap_mass_co2)    / giga / 1e3);
    fprintf('Residual: %f Gtons\n', sum(btrap_mass_co2_res)  / giga / 1e3);
    fprintf('Dissolved: %f Gtons\n', sum(btrap_mass_co2_dis) / giga / 1e3);

    % arrays of data for output:
    tot_trap_capa_sum(i,1)      = sum(tot_trap_capa) / giga / 1e3;
    strap_mass_co2_sum(i,1)     = sum(strap_mass_co2) / giga / 1e3;
    btrap_mass_co2_res_sum(i,1) = sum(btrap_mass_co2_res)  / giga / 1e3;
    btrap_mass_co2_dis_sum(i,1) = sum(btrap_mass_co2_dis) / giga / 1e3;

end

output.tot_trap_capa_sum        = tot_trap_capa_sum;
output.strap_mass_co2_sum       = strap_mass_co2_sum;
output.strap_mass_co2           = strap_mass_co2;
output.btrap_mass_co2_res_sum   = btrap_mass_co2_res_sum;
output.btrap_mass_co2_dis_sum   = btrap_mass_co2_dis_sum;

    
    
end




