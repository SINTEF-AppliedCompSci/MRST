function displayUnits = DisplayUnits(model)
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
end