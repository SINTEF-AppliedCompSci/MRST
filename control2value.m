
function val = control2value(u,parameters)
% Convert parameter param in model to control vector
np = numel(parameters);



reel = 1;
for k = 1:np
    switch parameters{k}.distribution
        case 'cell' %parameter disstribution per cell
            switch parameters{k}.name
                case 'porevolume'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   Indx = parameters{k}.Indx;
                   m = length(Indx);
                   val{k,1}(:,1) =  u(reel:reel+m-1)...
                                      *(umax-umin)+umin;
                   reel = m + 1;
                case 'initSw'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   Indx = parameters{k}.Indx;
                   m = length(Indx);
                   val{k,1}(:,1) =  u(reel:reel+m-1)...
                                      *(umax-umin)+umin;
                   reel = m + 1;
                %case {'swl', 'swcr', 'swu', 'sowcr'}  
                otherwise
                   warning('Parameter %s is not implemented',param{k})
            end
        case  'connection'
            switch parameters{k}.name
                case 'transmissibility'
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       val{k,1}(i,1) =  u(reel + i-1)...
                                      *(umax-umin)+umin;                                  
                   end
                   reel = reel + numel(parameters{k}.Indx) ;                   
                case 'porevolume'
                    % Asing a porevolume to each chanel
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       val{k,1}(i,1) =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                case 'permeability'
                    % Asing a porevolume to each chanel
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       val{k,1}(i,1) =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                otherwise
                   warning('Parameter %s is not implemented',param{k})
            end   
       case  'general'  
            switch parameters{k}.name
                case 'conntrans'
                   for i =  1 : size(parameters{k}.Indx,1)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       val{k,1}(i,1) =  u(reel + i-1)...
                                      *(umax-umin)+umin;      
                   end
                   reel = reel + size(parameters{k}.Indx,1) ;
                   
                case 'permeability'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   I_perm=  parameters{k}.Indx;
                       val{k,1}(1) =  u(reel )...
                                      *(umax-umin)+umin; 
                   reel = reel + 1 ;                                  
                otherwise
                   warning('Parameter %s is not implemented',param{k})
            end   
      otherwise
           warning('Parameter distribution %s is not implemented',paramDist{k})
      end
end
%val=val';
% Concatenate in a column vector

