function y = interpRegular2DTable(T, x, y)
    in_range_x = ( min(T.x)<=min(double(x)) && max(double(x)) <= max(T.x));
    in_range_y = ( min(T.y)<=min(double(y)) && max(double(y)) <= max(T.y));
    if(~(in_range_x && in_range_y))
       warning('value out of range of tensor table') 
    end
    
    if(isa(x,'ADI') || isa(y,'ADI'))
        xi=double(x);vi=double(y);
        [yi, dyidxi, dyidvi] = interp2DTableLoc(T, xi, vi);
        tmp={};
        if(isa(x,'ADI'))
         tmp{end+1} = lMultdiag(dyidxi, x.jac);
        end
        if(isa(y,'ADI'))
         tmp{end+1} = lMultdiag(dyidvi, y.jac);
        end
        if(numel(tmp)==1)
           y = ADI(yi, tmp{1});
           return
        elseif( numel(tmp)==2)   
            for i=1:numel(tmp{2})
                tmp{1}{i}=tmp{1}{i}+tmp{2}{i};
            end
        else
           error('Not possible')
        end                       
        y = ADI(yi, tmp{1});
        
    else
       y = interp2DTableLoc(T, x, y);
    end
end
function [yi, dyidxi, dyidvi] = interp2DTableLoc(table, P, T)
%function res = extract_val_vectorized(grid, spanP, spanT, P, T)
% Vectorized version of 'extract_val', where P and T are 
% allowed to be vectors
    [p_ix, t_ix, p, t, nans] = ix_and_local_par(table, P, T);
    
    
    lines = size(table.data, 1);
    corners = [table.data(p_ix   + lines * (t_ix - 1)), ...
               table.data(p_ix+1 + lines * (t_ix - 1)), ...
               table.data(p_ix   + lines * (t_ix    )), ...
               table.data(p_ix+1 + lines * (t_ix    ))];
  
    weight = [(1-p).*(1-t), p.*(1-t), (1-p).*t, p.*t];
                   
    yi = nan(numel(P), 1);
    yi(~nans) = sum(corners .* weight, 2);
    if(nargout>1)
        %dSpanP = diff(table.spanX);
        %dSpanT = diff(table.spanY);
        %dx = dSpanP/table.stepsX;
        %dy = dSpanT/table.stepsY;
        dx=table.dx;
        dy=table.dy;
        wx=[-(1-t),(1-t),-t,t];
        wy=[-(1-p),-p,(1-p),p];
        dyidxi=nan(numel(P), 1);
        dyidvi=nan(numel(P), 1);
        dyidxi(~nans)=sum(corners .* wx, 2)/dx;
        dyidvi(~nans)=sum(corners .* wy, 2)/dy;
    end
    
end
function [p_ix, t_ix, p, t, nans] = ix_and_local_par(table, P, T)
    
    nans = isnan(P) | isnan(T);
    
    Psteps = size(table.data, 1) - 1;
    Tsteps = size(table.data, 2) - 1;
    
    
    Dp=table.dx;
    Dt=table.dy;
    dSpanP=Dp*(Psteps);
    dSpanT=Dt*(Tsteps);
    p_ix = 1 + floor(((P(~nans) - table.x(1)) * Psteps)/dSpanP);    
    t_ix = 1 + floor(((T(~nans) - table.y(1)) * Tsteps)/dSpanT);
    
    p_ix = min(p_ix, Psteps);
    t_ix = min(t_ix, Tsteps);
    
    Psampled = table.x(1) + (p_ix - 1) * Dp;
    Tsampled = table.y(1) + (t_ix - 1) * Dt;
    
    p = (P(~nans) - Psampled) / Dp; % should be between 0 and 1
    t = (T(~nans) - Tsampled) / Dt; % should be between 0 and 1
    
    %assert(all(t<=1 & t>=0))
    %assert(all(p<=1 & p>=0))
end


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
