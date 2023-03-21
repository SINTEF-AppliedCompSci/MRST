function [h, h_max, dh] = saturation2Height(sg, sgmax, Gt, rg, rw, varargin)
%
% Considering a sharp interface system, compute the plume thickness (h) and
% historically maximum thickness (h_max) from upscaled saturation (sg) and
% historically maximum upscaled saturation (sgmax).  If porosity is passed
% along as a parameter, the vertical heterogeneity of porosity will be taken
% into account.
% 
% SYNOPSIS:
%   function [h, h_max, dh] = saturation2Height(sg, sgmax, Gt, rg, rw, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   sg       - upscaled saturation (current)
%   sgmax    - upscaled saturation (historical maximum)
%   Gt       - top surface grid
%   rg       - residual gas saturation value
%   rw       - residual water saturation value
%   varargin - the following optional parameters can be passed along:
%              * poro - the fine-scale (3D) porosity field.  If passed along,
%                       vertical heterogeneities in porosity will be taken
%                       into account.
%              * pvol_columns - total pore volume of each column.  The total
%                               pore volume is needed whenver vertical
%                               heterogeneities are taken into account.  If
%                               not provided, it will be computed internally
%                               whenever 'poro' is provided.
%              * tol - tolerance to use when converting from saturation to
%                      height when taking vertical heterogeneities into account.
%
% RETURNS:
%   h     - plume thickness
%   h_max - historically maximum plume thickness
%   dh    - value of derivative of h with respect to sg

%{
Copyright 2009-2023 SINTEF Digital, Mathematics & Cybernetics.

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
    
    opt = merge_options(struct('poro'         , []    , ...
                               'pvol_columns' , []    , ...
                               'tol'          , 1e-5) , ...
                        varargin{:});
    
    H = Gt.cells.H;
    ismax = (value(sg) == value(sgmax));
    
    if isempty(opt.poro)
        % porosity assumed uniform, we can use the simple transformation
        % formula:
        % s * H = h * (1 - rw) + (h_max - h) * rg
        % s_max * H = h_max * (1 - rw)


        h_max = sgmax .* H ./ (1 - rw);
        h = (H .* sg - h_max .* rg) ./ (1 - rg - rw);
    
        % correct possible values of h out of range (should only be due to numerical
        % inaccuracies/roundoff)
        h = min(max(h, 0), H); 

        dh = H ./ (1 - rg - rw);

        % The derivative has a kink at this point; the formula below
        % gives the derivative in the positive direction.  @@ Can
        % this lead to complications in the solver?
        dh(ismax) = H(ismax) ./ (1 - rw);

    else
        % The relationship between upscaled saturation and height depend on
        % the vertical porosity profile.  If we note the total pore volume of
        % a column that is covered by a plume with thickness h as V(h), we
        % have that s = V(h) / Vtot  (where Vtot is the total pore volume of
        % the column).

        column_cellnums = diff(Gt.cells.columnPos);
        map_columns = Gt.columns.cells;
        poro_areas(map_columns) = opt.poro(map_columns) .* ...
            rldecode(Gt.cells.volumes, column_cellnums);
        poro_areas = poro_areas(:);
        
        Vh = @(h, rw) pvol(value(h), rw, poro_areas, Gt);

        if isempty(opt.pvol_columns)
            Vtot = Vh(H, 0);
        else
            Vtot = opt.pvol_columns; % if we have precomputed values
                                     % available, use them instead
        end

        [h_max, dhmdv] = Vinv(Vh, value(sgmax) .* Vtot, Vtot, rw, H, opt.tol);
        
        vhmaxres = pvol(h_max, 1-rg, poro_areas, Gt);
        
        [h, dhdv] = Vinv(Vh, max(value(sg) .* Vtot - vhmaxres, 0), ...
                         Vtot, rw + rg, H, opt.tol);

        % convert dh/dv to dh/dsg
        dh = dhdv .* Vtot;

        % ensure h is ADI if required
        if isa(sg, 'ADI')
            h = ADI(h, lMultDiag(dh, sg.jac));
        end
        if isa(sgmax, 'ADI')
            h_max = ADI(h_max, lMultDiag(dhdv .* Vtot, sgmax.jac));
        end
    end
            
end

function [v, dvdh] = pvol(h, rw, poro_areas, Gt)
    [v, dvdh] = integrateVertically(poro_areas, h, Gt);
    
    % scaling results with fraction of porevolume available for CO2
    v    = v    * (1-rw);
    dvdh = dvdh * (1-rw);
    
end

function [h, dh] = Vinv(vfun, vtarget, vmax, rw, H, tol)

    h = vtarget ./ vmax .* H ./ (1-rw); % initial guess
    
    [v, dvdh] = vfun(h, rw);
    
    res = v - vtarget;
    
    dh = 1 ./ dvdh; % dh/dv
    
    while max(abs(res)) > tol
        h = h - res .* dh;
        h = max(min(h, H), 0);
        [v, dvdh] = vfun(h, rw);
        
        res = v - vtarget;
        dh = 1./ dvdh;
    end
end

function J = lMultDiag(d, J1)
   n = numel(d); 
   D = sparse((1:n)', (1:n)', d, n, n); 
   J = cell(1, numel(J1)); 
   for k = 1:numel(J)
      J{k} = D * J1{k}; 
   end
end
