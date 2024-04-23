function fluid = addVERelpermSharpInterface(fluid, Gt, rock, varargin)
    % If 'dh' is nonzero, caprock rugosity is modeled using the
    % accretion layer model.   
    % type can be 'simple', or 'integrated'

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

    opt = struct('krw', [], ...
                 'krg', [], ...
                 'dh', 0, ... % scalar or field
                 'type', 'simple');

    opt = merge_options(opt, varargin{:});
    
    [krw, krg] = deal(ifelse(isempty(opt.krw), 1-fluid.res_gas, opt.krw), ...
                      ifelse(isempty(opt.krg), 1-fluid.res_water, opt.krg));
    
    [poro3D, perm3D] = deal([]);
    [perm_columns, pvol_columns] = deal([]);
    if strcmpi(opt.type, 'integrated')
        assert(isfield(rock, 'parent')); % need fine-scale properties from parent
        assert(isfield(Gt, 'parent'));% we need access to fine-scale cell volumes
        
        poro3D = rock.parent.poro;
        perm3D = rock.parent.perm;
        pvol_columns = accumarray(rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos)), ...
                                  Gt.parent.cells.volumes(Gt.columns.cells) .* poro3D(Gt.columns.cells));
        perm_H = integrateVertically(perm3D, Gt.cells.H, Gt);
    end
    
    % helper function to convert from saturation to plume thickness
    h_hmax = @(sg, varargin) compute_heights(sg, get_sGmax(sg, varargin{:}), Gt, fluid, poro3D, pvol_columns);
    
                                                
    % accretion layer (upper zone where horizontal permeability is considered zero, 
    % constituting a very simple model of the impact of caprock rugosity)
    krH_lost = accretion_layer_perm(min(opt.dh, Gt.cells.H), Gt, perm3D);
        
    % Add permeability according to chosen model
    if strcmpi(opt.type, 'simple')
        fluid.krG = @(sg, p, varargin) krG_simple(h_hmax(sg, varargin{:}), Gt, fluid, rock, krg, krH_lost);
        fluid.krW = @(sw, p, varargin) krW_simple(h_hmax(1-sw, varargin{:}), Gt, fluid, rock, krw, krH_lost);
    else
        fluid.krG = @(sg, p, varargin) ...
            krG_integrated(h_hmax(sg, varargin{:}), Gt, fluid, rock, perm3D, krg, krH_lost, perm_H);
        fluid.krW = @(sw, p, varargin) ...
            krW_integrated(h_hmax(1-sw, varargin{:}), Gt, fluid, rock, perm3D, krw, krH_lost, perm_H);
    
    end        
    
    % Add capillary pressure
    fluid.pcWG = @(sg, p, varargin) pcWG(sg, p, fluid, Gt, poro3D, pvol_columns, varargin{:});
    fluid.invPc3D = @(p) invPc3D(p, fluid.res_water);    
end

% ----------------------------------------------------------------------------
function kr = krG_integrated(heights, Gt, fluid, rock, perm3D, kr_endpoint, krH_lost, perm_H)
    [H, h, hmax] = deal(Gt.cells.H, heights{:});
    assert(all(h >= 0));
    if isa(h, 'ADI')
        [krH_tmp, dkrH_tmp] = integrateVertically(perm3D(:,1), h.val, Gt);
        krH = ADI(krH_tmp, lMultDiag(dkrH_tmp, h.jac));
    else
        krH = integrateVertically(perm3D(:,1), h, Gt);
    end
    assert(all(krH >= 0));
    
    % deduct the impact of the accretion layer
    krH = max(krH - krH_lost, 0);
    
    % compute fraction of full permeability, and multiply by endpoint scaling.
    kr = krH ./ perm_H .* kr_endpoint;
        
end

% ----------------------------------------------------------------------------
function kr = krW_integrated(heights, Gt, fluid, rock, perm3D, kr_endpoint, krH_lost, perm_H)
    [H, h, hmax] = deal(Gt.cells.H, heights{:});
    assert(all(h >= 0));
    
    % computing permeability integrated along CO2 column (no modification for endpoint or 
    % caprock rugosity
    kr_H_co2 = krG_integrated(heights, Gt, fluid, rock, perm3D, 1.0, 0.0, perm_H) .* perm_H;
    
    kr_H_wat = perm_H - max(kr_H_co2, krH_lost);
    
    kr = kr_H_wat ./ perm_H .* kr_endpoint;
end

% ----------------------------------------------------------------------------
function krH = accretion_layer_perm(dh, Gt, perm3D)

    if isempty(perm3D)
        % simple, nonintegrated model
        krH = dh;
    else
        if isa(dh, 'ADI')
            [krH, dkrH] = integrateVertically(perm3D(:,1), dh.val, Gt);
            krH = ADI(krH, lMultDiag(dkrH, dh.jac));
        else
            krH = integrateVertically(perm3D(:,1), dh, Gt);
        end
    end
   assert(all(krH >= 0)); 
end

% ----------------------------------------------------------------------------
function res = ifelse(cond, yes, no)
    if cond
        res = yes;
    else
        res = no;
    end
end

% ----------------------------------------------------------------------------
function kr = krG_simple(heights, Gt, fluid, rock, kr_endpoint, krH_lost)
    [H, h, hmax] = deal(Gt.cells.H, heights{:});
    
    kr = max(h - krH_lost, 0) .* kr_endpoint ./ H;
end

% ----------------------------------------------------------------------------
function kr = krW_simple(heights, Gt, fluid, rock, kr_endpoint, krH_lost)
    [H, h, hmax] = deal(Gt.cells.H, heights{:});
    
    kr = (H - max(hmax, krH_lost) + ... % contribution to permeability, fully water-saturated zone
         max((hmax - max(h, krH_lost)), 0) * kr_endpoint) ./ ... % zone with residual co2 saturation
         H; % divide by full thickness
end

% ----------------------------------------------------------------------------
function heights = compute_heights(sg, sGmax, Gt, fluid, poro3D, pvol_columns)
    [h, hmax] = upscaledSat2height(sg, sGmax, Gt, ...
                                   'resSat', [fluid.res_water, fluid.res_gas], ...
                                   'poro', poro3D, ...
                                   'pvol_columns', pvol_columns);
    assert(all(h >= 0));
    heights = {h, hmax};
end

% ----------------------------------------------------------------------------
function sGmax = get_sGmax(sg, varargin)
    [opt, extra] = merge_options(struct('sGmax', []), varargin{:}); 
    sGmax = ifelse(isempty(opt.sGmax), sg, opt.sGmax);
end

% ----------------------------------------------------------------------------
function T = get_temperature(varargin)
    [opt, extra] = merge_options(struct('T', []), varargin{:});
    T = opt.T; % empty array if not provided
end


% ----------------------------------------------------------------------------
function s = invPc3D(p, res_water)
% return water saturation that is either 1 (if below sharp interface) or
% res_water (if above it).  @@ Note that from capillary pressure alone, it is 
% impossible to determine whether there is a zone of residual CO2 present,
% so residual CO2 is ignored here.
    s = (sign(p + eps) + 1) / 2 * (1 - res_water);  % CO2 saturation
    s = 1 - s; % converting to water saturation
end

% ----------------------------------------------------------------------------
function pc = pcWG(sg, p, fluid, Gt, poro3D, pvol_columns, varargin)
    
    sGmax = get_sGmax(sg, varargin{:});
    T = get_temperature(varargin{:});
    
    [h, h_max] = upscaledSat2height(sg, sGmax, Gt, ...
                                    'resSat', [fluid.res_water, fluid.res_gas], ...
                                    'poro', poro3D, ...
                                    'pvol_columns', pvol_columns);
     assert(all(h >= 0));

     if isempty(T)
        drho = (fluid.rhoWS .* fluid.bW(p) - fluid.rhoGS .* fluid.bG(p));
     else
         drho = (fluid.rhoWS .* fluid.bW(p, T) - fluid.rhoGS .* fluid.bG(p, T));
     end
     
     pc = drho .* h * norm(gravity);
    
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
