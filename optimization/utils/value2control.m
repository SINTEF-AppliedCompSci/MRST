function u = value2control(val,parameters)

% Convert parameter param in model to control vector
np = numel(parameters);

for k = 1:np
    switch parameters{k}.distribution
        case 'cell' %parameter disstribution per cell
            switch parameters{k}.name
                case {'porevolume','initSw','transmissibility'}
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                  Indx = parameters{k}.Indx;
                  ui   =  (val{k}(:)-umin)./(umax-umin);  
                  u{k} = ui;  
                  ui   =[];
                otherwise
                         error('Parameter %s is not implemented',parameters{k}.name)
            end
       case  'connection' 
            switch parameters{k}.name
                case {'transmissibility','porevolume','permeability','conntrans'} 
                  for i =  1 : numel(parameters{k}.Indx)
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      Indx = parameters{k}.Indx{i};
                      ui(i,1)   = (val{k}(i)-umin)./(umax-umin);
                  end                  
                  u{k} = ui;  
                  ui   =[];                 
                otherwise
                         error('Parameter %s is not implemented',parameters{k}.name)
            end
       case 'general'
             switch parameters{k}.name
                case {'transmissibility','porevolume','permeability'}
                  [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                  ui   = val{k};
                  u{k} = (ui-umin)./(umax-umin); 
                  ui   =[];
               case 'conntrans'                    
                   for i =  1 : size(parameters{k}.Indx,1)    
                      [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                      ui(i,1)   = (val{k}(i)-umin)./(umax-umin);
                  end                  
                  u{k} = ui;  
                  ui   =[]; 
                otherwise
                         error('Parameter %s is not implemented',parameters{k}.name)
            end
       otherwise
            error('Parameter distribution %s is not implemented',parameters{k}.distribution)
    end
end
% Concatenate in a column vector
u = vertcat(u{:});