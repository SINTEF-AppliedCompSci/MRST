function model = setrhoGmultfun(model, varargin)
%Undocumented Utility Function

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

    val = varargin{end};

    % The following functions in the fluid object must change:
    % - rhoG
    % - bG

    % At first call, we modify the fluid object to accomodate for the
    % changes
    if ~isfield(model.fluid, 'rhoGmult')
        % Keeping reference values of all fields that will be impacted by 
        % changes to 'rhoGmult'
        model.fluid.rhoGmult = 1;
        model.fluid.rhoGS_orig = model.fluid.rhoGS;
        model.fluid.rhoG_orig = model.fluid.rhoG;
        model.fluid.pcWG_orig = model.fluid.pcWG;
        model.fluid.bG_orig = model.fluid.bG;
    end
    
    % set multiplier and update density functions
    model.fluid.rhoGmult = val;
    model.fluid.rhoGS = model.fluid.rhoGS_orig .* val;
    model.fluid.rhoG = @(p, varargin) val .* model.fluid.rhoG_orig(p, varargin{:});
    
    % change in upscaled capillary pressure is best expressed through an
    % intermediary function
    cap_adjust = @(p) (model.fluid.rhoW(p) - val .* model.fluid.rhoG_orig(p)) ./ ...
                      (model.fluid.rhoW(p) - model.fluid.rhoG_orig(p));
    
    model.fluid.pcWG = @(sg, p, varargin) cap_adjust(p) .* ...
                                          model.fluid.pcWG_orig(sg, p, varargin{:});
end
