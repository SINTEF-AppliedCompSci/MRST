
function [model,schedule,state0] = control2problem(u,model,schedule,state0, parameters)
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
                   model.operators.pv(Indx) = u(reel:reel+m-1)...
                                        *(umax-umin)+umin;
                   reel = reel + m ;
                case 'initSw'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   Indx = parameters{k}.Indx;
                   m = length(Indx);
                   state0.s(Indx,1) = u(reel:reel+m-1)...
                                        *(umax-umin)+umin;
                   state0.s(Indx,2) = (1-u(reel:reel+m-1))...
                                        *(umax-umin)+umin;                 
                   reel = reel + m ;                   
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
                       model.operators.T(Indx) =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                case 'porevolume'
                    % Asing a porevolume to each chanel
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       model.operators.pv(Indx) =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                case 'permeability'
                    % Asing a porevolume to each chanel
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       if (parameters{k}.log == true)
                            model.rock.perm(Indx,1) = 10.^(u(reel + i-1)...
                                                          *(umax-umin)+umin);
                       else 
                            model.rock.perm(Indx,1) =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                       end
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                otherwise
                   warning('Parameter %s is not implemented',param{k})
            end
      case  'general'    
             switch parameters{k}.name
                case 'transmissibility'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   I_Tr =  parameters{k}.Indx;
                   model.operators.T(I_Tr) = 0*model.operators.T(I_Tr) + u(k)...
                                       *(umax-umin)+umin;  
                 case 'porevolume'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   I_pv =  parameters{k}.Indx;
                   model.operators.pv(I_pv) = 0*model.operators.pv(I_pv) + u(k)...
                                       *(umax-umin)+umin;
                 case 'permeability'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   I_perm=  parameters{k}.Indx;
                   model.rock.perm = 0*model.rock.perm + u(reel)...
                                       *(umax-umin)+umin; 
                   reel = reel + 1;
                case 'conntrans'
                   Indx = parameters{k}.Indx;
                   for i = 1 : size(parameters{k}.Indx,1)                      
                   [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                    for kk = 1 : numel( schedule.control)
                       schedule.control(kk).W(Indx(i,1)).WI(Indx(i,2)) = u(reel + i-1)...
                                                         *(umax-umin)+umin;
                    end
                   end
                   reel = reel + size(parameters{k}.Indx,1);
                 otherwise
                   warning('Parameter %s is not implemented',param{k}) 
             end
      otherwise
           warning('Parameter distribution %s is not implemented',paramDist{k})
      end
end
% Concatenate in a column vector

