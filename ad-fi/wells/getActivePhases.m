function actPh = getActivePhases(model)

    switch model
        
      case 'OW'
        actPh = [1 2];
      case 'WG'
        actPh = [1 3];
      case {'3P', 'BO', 'VO'}
        actPh = [1 2 3];
      otherwise
        error('Model not supported.');
    
    end

end