function displayUnits = DisplayUnits(model)
%
% DESCRIPTION: reads the plotting and output units from the model
%
% SYNOPSIS:
%   displayUnits = DisplayUnits(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - plot: plotting units
%
% RETURNS:
%   displayUnits - structure with the units of the plotting output in
%   string
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
displaySat    = 'fraction';
if(isfield(model.plot,'displaySat'))
    displaySat = model.plot.displaySat.inputUnit;
end
displayUnits.displaySat = displaySat;

displayTime   = 's';
if(isfield(model.plot,'displayTime'))
    displayTime = model.plot.displayTime.inputUnit;
end
displayUnits.displayTime = displayTime;

displayLength = 'm';
if(isfield(model.plot,'displayLength'))
    displayLength = model.plot.displayLength.inputUnit;
end
displayUnits.displayLength = displayLength;

displayVolume = 'm^3';
if(isfield(model.plot,'displayVolume'))
    displayVolume = model.plot.displayVolume.inputUnit;
end
displayUnits.displayVolume = displayVolume;

displayPress  = 'Pa';
if(isfield(model.plot,'displayPress'))
    displayPress = model.plot.displayPress.inputUnit;
end
displayUnits.displayPress = displayPress;

displayRate  = 'm^3/s';
if(isfield(model.plot,'displayRate'))
    displayRate = model.plot.displayRate.inputUnit;
end
displayUnits.displayRate = displayRate;