function appUnit = AppUnit(unit)
% ask siroos for the docs
    unit = lower(unit);
    
    % time    
    if(strcmp(unit,strcat('second')) || strcmp(unit,strcat('sec')) || strcmp(unit,strcat('s'))) 
        appUnit = 's'; 
        return;
    end    
    if(strcmp(unit,strcat('minute')) || strcmp(unit,strcat('min')) || strcmp(unit,strcat('m'))) 
        appUnit = 'min'; 
        return;
    end    
    if(strcmp(unit,strcat('hour')) || strcmp(unit,strcat('hr'))  || strcmp(unit,strcat('h'))) 
        appUnit = 'hr'; 
        return;
    end
    if(strcmp(unit,strcat('days')) || strcmp(unit,strcat('day'))  || strcmp(unit,strcat('d'))) 
        appUnit = 'day'; 
        return;
    end
    
    % length
    if(strcmp(unit,strcat('meter')) || strcmp(unit,strcat('m'))) 
        appUnit = 'm'; 
        return;
    end
    if(strcmp(unit,strcat('centimeter')) || strcmp(unit,strcat('cm'))) 
        appUnit = 'cm'; 
        return;
    end
    
    % volume
    if(strcmp(unit,strcat('m^3')) || strcmp(unit,strcat('m3'))) 
        appUnit = 'm^3'; 
        return;
    end
    if(strcmp(unit,strcat('cm^3')) || strcmp(unit,strcat('cm3'))) 
        appUnit = 'cm^3'; 
        return;
    end
    
    % volume rate
    if(strcmp(unit,strcat('m^3/second')) || strcmp(unit,strcat('m3/second')) || ...
       strcmp(unit,strcat('m^3/sec'))    || strcmp(unit,strcat('m3/sec'))    || ...
       strcmp(unit,strcat('m^3/s'))    || strcmp(unit,strcat('m3/s')))
        appUnit = 'm^3/s'; 
        return;
    end
    if(strcmp(unit,strcat('m^3/minute')) || strcmp(unit,strcat('m3/minute')) || ...
       strcmp(unit,strcat('m^3/min'))    || strcmp(unit,strcat('m3/min')))
        appUnit = 'm^3/min'; 
        return;
    end
    if(strcmp(unit,strcat('m^3/hour')) || strcmp(unit,strcat('m3/hour')) || ...
       strcmp(unit,strcat('m^3/hr'))   || strcmp(unit,strcat('m3/hr'))   || ...
       strcmp(unit,strcat('m^3/h'))    || strcmp(unit,strcat('m3/h')))
        appUnit = 'm^3/hr'; 
        return;
    end
    
    if(strcmp(unit,strcat('cm^3/second')) || strcmp(unit,strcat('cm3/second')) || ...
       strcmp(unit,strcat('cm^3/sec'))    || strcmp(unit,strcat('cm3/sec'))    || ...
       strcmp(unit,strcat('cm^3/s'))    || strcmp(unit,strcat('cm3/s')))
        appUnit = 'cm^3/s'; 
        return;
    end
    if(strcmp(unit,strcat('cm^3/minute')) || strcmp(unit,strcat('cm3/minute')) || ...
       strcmp(unit,strcat('cm^3/min'))    || strcmp(unit,strcat('cm3/min')))
        appUnit = 'cm^3/min'; 
        return;
    end
    if(strcmp(unit,strcat('cm^3/hour')) || strcmp(unit,strcat('cm3/hour')) || ...
       strcmp(unit,strcat('cm^3/hr'))   || strcmp(unit,strcat('cm3/hr'))   || ...
       strcmp(unit,strcat('cm^3/h'))    || strcmp(unit,strcat('cm3/h')))
        appUnit = 'cm^3/hr'; 
        return;
    end
    
    % pressure
    if(strcmp(unit,strcat('pa')) || strcmp(unit,strcat('pascal')))
        appUnit = 'Pa'; 
        return;
    end
    if(strcmp(unit,strcat('bar')) || strcmp(unit,strcat('bars')))
        appUnit = 'bar'; 
        return;
    end
    if(strcmp(unit,strcat('atm')) || strcmp(unit,strcat('atmospher')))
        appUnit = 'atm'; 
        return;
    end
    
    error('Input unit could not be found, check the input unit.');
end