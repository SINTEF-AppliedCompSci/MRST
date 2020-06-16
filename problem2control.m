

function u = problem2control(model,schedule,state0,parameters)

% Convert parameter param in model to control vector
np = numel(parameters);




for k = 1:np
    switch parameters{k}.distribution
        case 'cell' %parameter disstribution per cell
            switch parameters{k}.name
                case 'porevolume'
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                  Indx = parameters{k}.Indx;
                  ui   = model.operators.pv(Indx);
                  u{k} = (ui-umin)./(umax-umin);  
                  ui   =[];
                case 'initSw'
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                  Indx = parameters{k}.Indx;
                  ui   = state0.s(Indx,1);
                  u{k} = (ui-umin)./(umax-umin);  
                  ui   =[];
                otherwise
                         warning('Parameter %s is not implemented',params{k})
            end
       case  'connection' 
            switch parameters{k}.name
                case 'transmissibility' 
                  for i =  1 : numel(parameters{k}.Indx)
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      Indx = parameters{k}.Indx{i};
                      ui(i,1)   = (mean(model.operators.T(Indx))-umin)./(umax-umin);
                  end                  
                  u{k} = ui;  
                  ui   =[];
                case 'porevolume'                 
                  for i =  1 : numel(parameters{k}.Indx)
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      Indx = parameters{k}.Indx{i};
                      ui(i,1)   = (mean(model.operators.pv(Indx))-umin)./(umax-umin);
                  end                  
                  u{k} = ui;   
                  ui   =[];
                case 'permeability'                 
                  for i =  1 : numel(parameters{k}.Indx)
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      Indx = parameters{k}.Indx{i};
                      if (parameters{k}.log == true)
                         ui(i,1)   = (log10(mean(model.rock.perm(Indx,1)))-umin)./(umax-umin);
                      else
                         ui(i,1)   = (mean(model.rock.perm(Indx,1))-umin)./(umax-umin);
                      end
                  end                  
                  u{k} = ui;   
                  ui   =[];
                otherwise
                         warning('Parameter %s is not implemented',params{k})
            end
       case 'general'
             switch parameters{k}.name
                case  {'swl', 'swcr', 'swu', 'sowcr'} 
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                 
                  ui   = model.operators.T(1);
                  u{k} = (ui-umin)./(umax-umin); 
                  ui   =[];
                case 'porevolume'
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   ui   = model.operators.pv(2);
                   u{k} = (ui-umin)./(umax-umin);
                   ui   =[];
                case 'permeability'
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   ui   = model.rock.perm(1);
                   u{k} = (ui-umin)./(umax-umin);
                   ui   =[];
                case 'conntrans'
                   Indx = parameters{k}.Indx;                   
                   for i =  1 : size(parameters{k}.Indx,1) 
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      ui(i,1)   = (schedule.control.W(Indx(i,1)).WI(Indx(i,2))-umin)./(umax-umin);
                   end
                   u{k} = ui;
                   ui   = [];
                  end 
                otherwise
                         warning('Parameter %s is not implemented',params{k})
            end
       otherwise
            warning('Parameter distribution %s is not implemented',paramDist{k})
    end
end
% Concatenate in a column vector
u = vertcat(u{:});
