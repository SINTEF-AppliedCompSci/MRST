function [model,schedule,state0] = control2problem(u,model,schedule,state0, parameters)
% Convert parameter param in model to control vector

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

np = numel(parameters);

reel = 1;
for k = 1:np
    switch parameters{k}.distribution
        case 'cell' %parameter disstribution per cell
            switch parameters{k}.name
                case 'transmissibility'
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   Indx = parameters{k}.Indx;
                   m = length(Indx);
                   model.operators.T(Indx) = u(reel:reel+m-1)...
                                        *(umax-umin)+umin;
                   reel = reel + m ;
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
                   error('Parameter %s is not implemented',param{k})
            end
        case  'connection'
            switch parameters{k}.name
                
                case {'transmissibility','porevolume'}
                    % handle mapping of names
                   switch parameters{k}.name
                       case 'transmissibility'
                          fname='T'
                       case 'porevolume'
                          fname='pv' 
                       otherwise
                           error()                           
                   end                    
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       val =  u(reel + i-1)*(umax-umin)+umin;
                       switch parameters{k}.type
                           case 'value'
                                model.operators.(fname)(Indx) =  val;
                           case 'multiplier'
                               model.operators.(fname)(Indx) = model.operators.(fname)(Indx)*val;
                           otherwise
                               error('Type not valid: %s ', parameters{k}.type);
                       end              
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                case 'conntrans'                  
                   for i = 1 : numel(parameters{k}.Indx)
                       Indx = parameters{k}.Indx{i};
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2));
                       val = u(reel + i-1)*(umax-umin)+umin;
                       for ii=1:size(Indx,1)
                          for kk = 1 : numel( schedule.control)
                               switch parameters{k}.type
                                    case 'value'
                                        schedule.control(kk).W(Indx(ii,1)).WI(Indx(ii,2)) = val;
                                    case 'multiplier'
                                        schedule.control(kk).W(Indx(ii,1)).WI(Indx(ii,2)) = schedule.control(kk).W(Indx(ii,1)).WI(Indx(ii,2))*val;
                                   otherwise
                                        error('Type not valid: %s ', parameters{k}.type);
                               end
                          end
                       end    
                   end                   
                   reel = reel + numel(parameters{k}.Indx);   
                case 'permeability'
                    % Asing a porevolume to each chanel
                   for i =  1 : numel(parameters{k}.Indx)
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2)); 
                       Indx = parameters{k}.Indx{i};
                       if (parameters{k}.log == true)
                            model.rock.perm(Indx,1) = 10.^(u(reel + i-1)...
                                                          *(umax-umin)+umin);
                       else 
                            val =  u(reel + i-1)...
                                      *(umax-umin)+umin;
                            switch parameters{k}.type
                                case 'value'
                                    model.rock.perm(Indx,1) =  val;
                                case 'multiplier'
                                    model.rock.perm(Indx,1) = model.rock.perm(Indx,1)*val;
                                otherwise
                                    error('Type not valid: %s ', parameters{k}.type);
                            end             
                       end                       
                   end
                   reel = reel + numel(parameters{k}.Indx) ;
                otherwise
                   error('Parameter %s is not implemented',param{k})
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
                 case {'swl','swcr', 'swu', 'sgl', ...
                       'sgcr','sgu','sowcr','sogcr',...
                       'krw','kro','krg'}
                   [umin, umax] = deal(parameters{k}.boxLims(1), parameters{k}.boxLims(2)); 
                   scaling = {upper(parameters{k}.name), u(k)...
                                                        *(umax-umin)+umin};
                   model = imposeRelpermScaling(model, scaling{:});
                case 'conntrans'
                   Indx = parameters{k}.Indx;
                   for i = 1 : size(parameters{k}.Indx,1)                      
                       [umin, umax] = deal(parameters{k}.boxLims(i,1), parameters{k}.boxLims(i,2));
                       val = u(reel + i-1)*(umax-umin)+umin;
                       for kk = 1 : numel( schedule.control)
                           switch parameters{k}.type
                            case 'value'
                                schedule.control(kk).W(Indx(i,1)).WI(Indx(i,2)) = val;
                            case 'multiplier'
                                schedule.control(kk).W(Indx(i,1)).WI(Indx(i,2)) = schedule.control(kk).W(Indx(i,1)).WI(Indx(i,2))*val;
                            otherwise
                               error('Type not valid: %s ', parameters{k}.type);
                           end
                        end
                   end                   
                   reel = reel + size(parameters{k}.Indx,1);
                 otherwise
                   error('Parameter %s is not implemented',param{k}) 
             end
      otherwise
           error('Parameter distribution %s is not implemented',parameters{k}.distribution)
      end
end
% Concatenate in a column vector
