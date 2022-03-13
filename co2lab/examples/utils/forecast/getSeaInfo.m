function [ info ] = getSeaInfo( name, rhoCref )
% Returns relevant info for a specified sea or formation in the Norwegian
% Contential Shelf.
%
% SYNOPSIS:
%   info = getSeaInfo('NorthSea');
%   info = getSeaInfo('Utsirafm');
%
% DESCRIPTION:
%   Information collected from NPD's Compiled CO2 Storage Atlas (2014), as
%   well as other sources such as the Sleipner benchmark, is returned for
%   one of the three Seas in the Norwegian Continental Shelf, or for a
%   particular formation name.
%
%   Density and viscosity does not refer to sea water properties, but
%   rather refers to the water located under the sea bottom in the
%   geological formation.
%
%   Units:
%       seafloor_depth  (meters)
%       seafloor_temp   (Celsius)
%       temp_gradient   (Celsius/kilometer)
%           NB: temp_gradient must be supplied in units of degrees C/km,
%           because unit conversion to meters is performed in local helper
%           functions.
%       rho             (kg/m3)
%       mu              (Pa*s)
%       res_sat         unitless
%       co2_solubility  (kg/m3)
%           NB: co2_solubility is a function of pressure and temperature,
%           thus could be computed specifically for each formation's
%           conditions.
%       dis_max         unitless
%       press_deviation (percent)

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

    northSeaNames       = getNorthSeaNames();
    norwegianSeaNames   = getNorwegianSeaNames();
    barentsSeaNames     = getBarentsSeaNames();
    

    if strcmpi(name,'NorthSea') || any(strcmpi(northSeaNames,name))
        
        % Atlas, chp 4.

        info.seafloor_depth  = 100 * meter; % @@ or should be 80 m for Utsira?
        info.seafloor_temp   =  7;
        info.temp_gradient   = 35.6;
        
        [rho, mu, sr, sw]    = getValuesSPE134891(); % from Sleipner benchmark
        info.water_density   = rho(1);
        info.co2_density     = rho(2);
        info.water_mu        = mu(1);
        info.co2_mu          = mu(2);
        info.res_sat_co2     = sr; 
        info.res_sat_wat     = sw;
        
        info.co2_solubility  = 53 * kilogram / meter^3;       % Chadwick et al 2008 (for Utsira, in brine?)
        info.dis_max         = info.co2_solubility / rhoCref;
        
        info.press_deviation = 0; % pressure devation from hydrostatic (%)
        
        info.rhoCref = rhoCref;
        info.sea = 'North';

        
        
    elseif strcmpi(name,'NorwegianSea') || any(strcmpi(norwegianSeaNames,name))
        
        % Atlas, chp 5.
        %warning('Some sea values need to be confirmed.')

        info.seafloor_depth  = 225 * meter;  % varies between 100-500 m in depth, but generally 200-250 m (personal communication with Ane Lothe)
        info.seafloor_temp   =  5;           % Lundin et al 2005 report
        info.temp_gradient   = 41.3;         % Lundin et al 2005 report
        % NB: Lothe et al 2014 used 40 C/km in their simulations
        
        [rho, mu, sr, sw]    = getValuesSPE134891();
        info.water_density   = rho(1);      % NB: confirm @@
        info.co2_density     = 700;         % Atlas, chp 5, pg 105
        info.water_mu        = mu(1);       % NB: confirm @@
        info.co2_mu          = mu(2);       % NB: confirm @@
        info.res_sat_co2     = sr;          % NB: confirm @@
        info.res_sat_wat     = sw;          % NB: confirm @@
        
        info.co2_solubility  = 53 * kilogram / meter^3;       % Chadwick et al 2008 (for Utsira, in brine?)
        info.dis_max         = info.co2_solubility / rhoCref;
        
        info.press_deviation = 0;
        
        info.rhoCref = rhoCref;
        info.sea = 'Norwegian';
        
        
    elseif strcmpi(name,'BarentsSea') || any(strcmpi(barentsSeaNames,name))
        
        % Atlas, chp 6.
        %warning('Some sea values need to be confirmed.')
        
        % NB: chp 6, pg 128 in Atlas says water density is 1.1 g/cm3
        % and described it as strongly saline (> 100 000 ppm) in the
        % Snohvit Field.
        
        % NB: CO2 density in Bjarmeland Platform is reportedly 650 kg/m3,
        % and 700 kg/m3 in other formations.
        %warning('CO2 density in Bjarmeland Platform is reportedly 650 kg/m3.')
        
        % Info on seafloor_depth and seafloor_temp:
        %   Atlas, chp 2 states Barents sea is up to 500 meter deep, with
        %   sea floor temperatures of 0 degrees Celsius or lower. However,
        %   the Atlas, chp 6, pg 134 reports that the sea depth at the
        %   Snohvit field is 330 meters.
        %
        %   See map and reference in Atlas, chp 6, pg 112. Sea bathymetry
        %   appears to be around 400-600 meters in area of Hammerfest.
        %
        %   Using a temperature gradient of 40 C/km (from temp gradient map
        %   in Atlas, chp 3, pg 27), a seafloor temp of 4 C was calculated
        %   to meet a constraint of temp=98 C at a depth of 2.6km (which
        %   was given in Pham et al 2011 for the Tubaen formation).
        %
        %   See also: temp gradient data in fig 6, Halland and Riis 2014.
        
        info.seafloor_depth  = 330 * meter;  % reported for Snohvit, Atlas, chp 6, pg 134
        info.seafloor_temp   =  4;           % computed to meet conditions reported in Pham et al 2011 for Tubaen
        info.temp_gradient   = 40;           % Atlas, chp 3, pg 27
        
        [~, mu, sr, sw]      = getValuesSPE134891();
        info.water_density   = 1100;        % Atlas, chp 6, pg 128
        info.co2_density     = 700;         % Atlas, chp 6, pg 144
        info.water_mu        = mu(1);       % confirm @@
        info.co2_mu          = mu(2);       % confirm @@
        info.res_sat_co2     = sr;
        info.res_sat_wat     = sw;
        
        info.co2_solubility  = 53 * kilogram / meter^3;       % Chadwick et al 2008 (for Utsira, in brine?)
        info.dis_max         = info.co2_solubility / rhoCref;
        
        info.press_deviation = 0;
        
        info.rhoCref = rhoCref;
        info.sea = 'Barents';
        
        
    else
        error(['Unknown sea or formation name. ', ...
            'Sea name options are: NorthSea, NorwegianSea, BarentsSea. ', ...
            'Formation name options are: '])
    end
    
end
