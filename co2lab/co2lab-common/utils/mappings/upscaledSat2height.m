function [h, h_max, dh] = upscaledSat2height(S, S_max, Gt, varargin)
% Compute upscaled height, based on upscaled saturation
%
% SYNOPSIS:
%   function [h, h_max] = upscaledSat2height(S, S_max, Gt, varargin)
%
% DESCRIPTION:
% Compute upscaled height (present and max), based on upscaled saturation
% (present and max), based either on a sharp-interface model, or a general
% model where the upscaled capillary pressure function is provided by the
% caller.
% 
% In the case of a sharp interface model, the derivative dh/ds can also be
% returned as an optional third output argument.
%    
% The sharp interface model also supports vertical heterogeneities in the
% porosity field.  In this case, the porosity field must be provided as an
% additional argument through `varargin` (see below).
%
% PARAMETERS:
%   S        - Upscaled present CO2 saturation
%   S_max    - Upscaled, historically maximum CO2 saturation
%   Gt       - Top-surface grid in question
%   varargin - Additional parameters (key, value), depending of conversion model:
%              * If a _sharp interface model_  is assumed, then the function 
%                needs the following argument:
%                - 'resSat' - [rw, rc], where 'rw' is the residual water
%                             saturation (assumed constant), and 'rc' is the
%                             residual CO2 saturation.
%                Moreover, if the vertical porosity field is not uniform, 
%                then the following additional arguments are needed:
%               - 'poro'         - Vertical porosity field (fine-scale, 3D)
%               - 'pvol_columns' - total pore volume of each column.  The total
%                                  pore volume is needed whenever vertical
%                                  heterogeneities are taken into account.  If
%                                  not provided, it will be computed internally
%                                  whenever 'poro' is provided.
%               - 'tol'          - tolerance to use when converting from 
%                                  saturation to height when taking vertical 
%                                  heterogeneities into account.
%   
%              * If a _general_ model is assumed, then the function needs the
%                following argument:
%                - pcWG(S, p, S_max) - upscaled capillary pressure as a
%                                      function of upscaled saturation,
%                                      current pressure and max. upscaled
%                                      saturation.  This function will be
%                                      available from the fluid object if a
%                                      capillary fringe model is used.
%                - rhoW(p)           - density of water [oil] phase, as a
%                                      function of pressure.
%                - rhoG(p)           - density of CO2 [gas] phase, as a
%                                      function of pressure.
%                - 'p'               - current pressure
%                  (Upscaled capillary pressure here defined as the
%                  pressure difference between phases at the level
%                  of the caprock, assuming that the difference is 0
%                  at depth 'h' (the deepest point where there is
%                  free flow of CO2).
%
% RETURNS:
%   h     - The height, corresponding to the vertical distance between
%           caprock and the deepest point for which there is still nonzero
%           CO2 flow.
%   h_max - The historically maximum height

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('resSat', [0, 0], 'pcWG', [], 'p', [], 'rhoW', [], 'rhoG', [], ...
                 'poro', [], 'pvol_columns', [], 'tol', 1e-5);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.pcWG)
        % assume a sharp interface model (possibly with vertical
        % heterogeneities)
        [rw, rg] = deal(opt.resSat(1), opt.resSat(2));
        [h, h_max, dh] = ...
            sharp_interface_case(S, S_max, Gt, rg, rw, opt.poro, ...
                                 opt.pvol_columns, opt.tol);
    else 
        % Capillary pressure function provided - we assume a general models
        assert(nargout < 3) % dh not currently supported for general model
        [h, h_max] = capillary_case(S, S_max, Gt, opt.pcWG, opt.p, ...
                                    opt.rhoW, opt.rhoG);
    end
end

% ----------------------------------------------------------------------------
function [h, h_max, dh] = sharp_interface_case(sg, sgmax, Gt, rg, rw, ...
                                               poro, pvol_columns, tol)
    
    H = Gt.cells.H;
    ismax = (value(sg) == value(sgmax));
    
    if isempty(poro)
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
        poro_areas(map_columns) = poro(map_columns) .* ...
            rldecode(Gt.cells.volumes, column_cellnums);
        poro_areas = poro_areas(:);
        
        Vh = @(h, rw) pvol(value(h), rw, poro_areas, Gt);

        if isempty(pvol_columns)
            Vtot = Vh(H, 0);
        else
            Vtot = pvol_columns; % if we have precomputed values
                                 % available, use them instead
        end

        [h_max, dhmdv] = Vinv(Vh, value(sgmax) .* Vtot, Vtot, rw, H, tol);
        
        vhmaxres = pvol(h_max, 1-rg, poro_areas, Gt);
        
        [h, dhdv] = Vinv(Vh, max(value(sg) .* Vtot - vhmaxres, 0), ...
                         Vtot, rw + rg, H, tol);

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

% ----------------------------------------------------------------------------
function [h, h_max] = capillary_case(S, S_max, Gt, pcWG, p, rhoW, rhoG)
    assert(~isempty(p));
    assert(~isempty(pcWG));
    assert(~isempty(rhoW));
    assert(~isempty(rhoG));
    pc    = pcWG(S, p, 'sGmax', S_max);
    pcmax = pcWG(S_max, p, 'sGmax', S_max);
    drho  = norm(gravity) * (rhoW(p) - rhoG(p));
    h     = pc ./ drho;
    h_max = pcmax ./ drho;    
end

% ----------------------------------------------------------------------------
function [v, dvdh] = pvol(h, rw, poro_areas, Gt)
    [v, dvdh] = integrateVertically(poro_areas, h, Gt);
    
    % scaling results with fraction of porevolume available for CO2
    v    = v    * (1-rw);
    dvdh = dvdh * (1-rw);
    
end

% ----------------------------------------------------------------------------
function [h, dh] = Vinv(vfun, vtarget, vmax, rw, H, tol)

    h = vtarget ./ vmax .* H ./ (1-rw); % initial guess
    
    [v, dvdh] = vfun(h, rw);
    
    res = v - vtarget;
    
    dh = 1 ./ dvdh; % dh/dv

    count = 0;
    while max(abs(res)) > tol
        count = count + 1;
        if count > 10
            error('Could not converge while computing h from volume.');
        end

        h = h - res .* dh;
        h = max(min(h, H), 0);
        [v, dvdh] = vfun(h, rw);
        
        res = v - vtarget;
        dh = 1./ dvdh;
    end
end

% ----------------------------------------------------------------------------
function J = lMultDiag(d, J1)
   n = numel(d); 
   D = sparse((1:n)', (1:n)', d, n, n); 
   J = cell(1, numel(J1)); 
   for k = 1:numel(J)
      J{k} = D * J1{k}; 
   end
end
