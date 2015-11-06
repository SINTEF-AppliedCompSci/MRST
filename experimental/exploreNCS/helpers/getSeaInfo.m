function [ info ] = getSeaInfo( seaName )
% Returns relevant info for a specified sea.
%
% SYNOPSIS:
%   info = getSeaInfo('NorthSea');
%
% DESCRIPTION:
%   Information collected from NPD's Compiled CO2 Storage Atlas (2011), as
%   well as other sources, is returned for one of the three Seas in the
%   Norwegian Continental Shelf.
%
%   Density and viscosity does not refer to sea body properties, but rather
%   refers to the geological formation water located sub-sea.
%
%

    if strcmpi(seaName,'NorthSea')
        
        % Atlas, chp 4.

        info.seafloor_depth = 100 * meter;
        info.seafloor_temp  =  7;           % in Celsius
        info.temp_gradient  = 35.6;         % degrees per kilometer
        
        [rho, mu, sr, sw] = getValuesSPE134891(); % from Sleipner benchmark
        info.water_density   = rho(1);      % kg per m3
        info.co2_density     = rho(2);
        info.water_mu        = mu(1);       % Pa per s
        info.co2_mu          = mu(2);
        info.res_sat_co2     = sr; 
        info.res_sat_wat     = sw;
        
        rhoCref      = info.co2_density;
        info.dis_max = (53 * kilogram / meter^3) / rhoCref; % from CO2store
        
        info.press_deviation = 0; % pressure devation from hydrostatic (%)

        % Optional:
        info.fmNames = getNorthSeaNames();
        
        
    elseif strcmpi(seaName,'NorwegianSea')
        
        % Atlas, chp 5.
        warning('Some sea values need to be confirmed.')

        info.seafloor_depth = 350 * meter;  % NB: confirm @@
        info.seafloor_temp  =  5;           % Lundin et al 2005 report
        info.temp_gradient  = 41.3;         % Lundin et al 2005 report
        
        [rho, mu, sr, sw] = getValuesSPE134891();
        info.water_density   = rho(1);      % NB: confirm @@
        info.co2_density     = 700;         % Atlas, chp 5, pg 105
        info.water_mu        = mu(1);       % NB: confirm @@
        info.co2_mu          = mu(2);       % NB: confirm @@
        info.res_sat_co2     = sr;          % NB: confirm @@
        info.res_sat_wat     = sw;          % NB: confirm @@
        
        rhoCref      = info.co2_density;
        info.dis_max = (53 * kilogram / meter^3) / rhoCref;
        
        info.press_deviation = 0;

        % Optional:
        info.fmNames = getNorwegianSeaNames();
        
        
    elseif strcmpi(seaName,'BarentsSea')
        
        % Atlas, chp 6.
        warning('Some sea values need to be confirmed.')
        
        % NB: chp 6, pg 128 in Atlas says water density is 1.1 g/cm3
        % and described it as strongly saline (> 100 000 ppm) in the
        % Snohvit Field.
        
        % NB: CO2 density in Bjarmeland Platform is reportedly 650 kg/m3,
        % and 700 kg/m3 in other formations.
        
        % Pham et al 2011 states Tubaen formation has an average temp of 98
        % C, at a depth of 2600 meters. A temp gradient of 40 C/km was
        % calculated to meet this constraint. See also: temp gradient data
        % in fig 6, Halland and Riis 2014
        
        info.seafloor_depth = 330 * meter;  % Atlas reference? @@
        info.seafloor_temp  =  4;           % Atlas reference? @@
        info.temp_gradient  = 40;
        
        [~, mu, sr, sw] = getValuesSPE134891();
        info.water_density   = 1100;        % Atlas, chp 6, pg 128
        info.co2_density     = 700;         % Atlas, chp 6, pg 144
        warning('CO2 density in Bjarmeland Platform is reportedly 650 kg/m3.')
        info.water_mu        = mu(1);
        info.co2_mu          = mu(2);
        info.res_sat_co2     = sr; 
        info.res_sat_wat     = sw;
        
        rhoCref      = info.co2_density;
        info.dis_max = (53 * kilogram / meter^3) / rhoCref;
        
        info.press_deviation = 0;
        
        % Optional:
        info.fmNames = getBarentsSeaNames();
        
        
    else
        error(['Unknown sea name. Options are: NorthSea, NorwegianSea,', ...
               ' BarentsSea.'])
    end
    
end

