function conFac = Convert(unit)
%
% DESCRIPTION: Reads the input unit as a string a returns the corresponding
%              conversion factor to SI units
%
% SYNOPSIS:
%   conFac = Convert(unit)
%
% PARAMETERS:
%   unit - string defining the unit
%
% RETURNS:
%   conFac - numerical value of conversion factor to SI units
%
% EXAMPLE:
%   Convert("cm");
%   Convert("rpm");
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
    unit = lower(unit);
    
    % volume rate
    splitUnit = split(string(unit),'/');
    if(length(splitUnit) > 1)
        conFac = Convert(splitUnit(1));
        for i = 2 : length(splitUnit)
            conFac = conFac / Convert(splitUnit(i));
        end
        return; 
    end
    
    % rpm, rotation
    if(strcmp(unit,strcat('rpm'))), conFac = 2 * pi() / 60; return; end
    
    % percent
    if(strcmp(unit,strcat('percent'))), conFac = centi; return; end
    if(strcmp(unit,strcat('%'))),       conFac = centi; return; end
    
    % fraction
    if(strcmp(unit,strcat('fraction'))), conFac = 1; return; end
    if(strcmp(unit,strcat('none'))),     conFac = 1; return; end
    if(strcmp(unit,strcat('-'))),        conFac = 1; return; end
    if(strcmp(unit,strcat(''))),         conFac = 1; return; end
    
    % time
    if(strcmp(unit,strcat('s'))),      conFac = second; return; end
    if(strcmp(unit,strcat('second'))), conFac = second; return; end
    if(strcmp(unit,strcat('s^2'))),    conFac = (second)^2; return; end
    if(strcmp(unit,strcat('s2'))),    conFac = (second)^2; return; end
    if(strcmp(unit,strcat('min'))),    conFac = minute; return; end
    if(strcmp(unit,strcat('minute'))), conFac = minute; return; end
    if(strcmp(unit,strcat('h'))),      conFac = hour; return; end
    if(strcmp(unit,strcat('hr'))),     conFac = hour; return; end
    if(strcmp(unit,strcat('hour'))),   conFac = hour; return; end
    if(strcmp(unit,strcat('d'))),      conFac = day; return; end
    if(strcmp(unit,strcat('day'))),    conFac = day; return; end
    if(strcmp(unit,strcat('days'))),   conFac = day; return; end
    
    % length
    if(strcmp(unit,strcat('m'))),  conFac = meter; return; end
    if(strcmp(unit,strcat('cm'))), conFac = centi*meter; return; end
    if(strcmp(unit,strcat('mm'))), conFac = milli*meter; return; end
    
    % volume
    if(strcmp(unit,strcat('m3'))),   conFac = (meter)^3; return; end
    if(strcmp(unit,strcat('m^3'))),  conFac = (meter)^3; return; end
    if(strcmp(unit,strcat('cm3'))),  conFac = (centi*meter)^3; return; end
    if(strcmp(unit,strcat('cm^3'))), conFac = (centi*meter)^3; return; end
    if(strcmp(unit,strcat('ft3'))),  conFac = (ft)^3; return; end
    if(strcmp(unit,strcat('ft^3'))), conFac = (ft)^3; return; end
    if(strcmpi(unit,strcat('ml'))), conFac = (milli*litre); return; end
    
%     volume rate ( fractions are not needed used anymore )
%     if(strcmp(unit,strcat('m3/s'))),       conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m^3/s'))),      conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m3/sec'))),       conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m^3/sec'))),      conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m3/second'))),  conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m^3/second'))), conFac = (meter)^3/second; return; end
%     if(strcmp(unit,strcat('m3/m'))),       conFac = (meter)^3/minute; return; end
%     if(strcmp(unit,strcat('m^3/m'))),      conFac = (meter)^3/minute; return; end
%     if(strcmp(unit,strcat('m3/min'))),  conFac = (meter)^3/minute; return; end
%     if(strcmp(unit,strcat('m^3/min'))), conFac = (meter)^3/minute; return; end
%     if(strcmp(unit,strcat('m3/minute'))),  conFac = (meter)^3/minute; return; end
%     if(strcmp(unit,strcat('m^3/minute'))), conFac = (meter)^3/minute; return; end

    
    % pressure
    if(strcmp(unit,strcat('pa'))),    conFac = Pascal; return; end
    if(strcmp(unit,strcat('psia'))),  conFac = psia; return; end
    if(strcmp(unit,strcat('bar'))),   conFac = barsa; return; end
    if(strcmp(unit,strcat('barsa'))), conFac = barsa; return; end
    if(strcmp(unit,strcat('atm'))),   conFac = atm; return; end
    
    % mass
    if(strcmp(unit,strcat('kg'))), conFac = kilogram; return; end
    if(strcmp(unit,strcat('g'))),  conFac = gram; return; end
    if(strcmp(unit,strcat('lb'))), conFac = pound; return; end
    
%     density ( fractions are not needed used anymore )
%     if(strcmp(unit,strcat('kg/m^3'))), conFac = kilogram/(meter)^3; return; end
%     if(strcmp(unit,strcat('kg/m3'))),  conFac = kilogram/(meter)^3; return; end
    
    % viscosity     
    if(strcmp(unit,strcat('cp'))),   conFac = centi*poise; return; end
    if(strcmp(unit,strcat('pa.s'))), conFac = Pascal*second; return; end
    
    % permeability
    if(strcmp(unit,strcat('md'))), conFac = milli*darcy; return; end
    if(strcmp(unit,strcat('d'))),  conFac = darcy; return; end
    
%     % acceleration
%     if(strcmp(unit,strcat('m/s^2'))), conFac = meter/(second)^2; return; end
    
    error('Input unit could not be found, check the input unit.');
end