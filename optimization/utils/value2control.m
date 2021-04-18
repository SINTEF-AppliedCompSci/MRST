function u = value2control(val,parameters)
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
