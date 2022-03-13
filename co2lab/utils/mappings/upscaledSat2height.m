function [h, h_max] = upscaledSat2height(S, S_max, Gt, varargin)
% Compute upscaled height, based on upscaled saturation
%
% SYNOPSIS:
%   function [h, h_max] = upscaledSat2height(S, S_max, varargin)
%
% DESCRIPTION:
% Compute upscaled height (present and max), based on upscaled saturation
% (present and max), based either on a sharp-interface model, or a general
% model where the upscaled capillary pressure function is provided by the
% caller.
%    
% PARAMETERS:
%   S        - Upscaled present saturation
%   S_max    - Upscaled, historically maximum, saturation
%   Gt       - Top-surface grid in question
%   varargin - Additional parameters (key, value), depending of conversion model:
%              * If a _sharp interface model_ (with vertically constant rock
%              properties) is assumed, then the function needs the following
%              argument:
%                - 'resSat' - [rw, rc], where 'rw' is the residual water
%                             saturation (assumed constant), and 'rc' is the
%                             residual CO2 saturation.
%              * If a _general_ model is assumed, then the function needs the
%              following argument:
%                - pcWG(S, p, S_max) - upscaled capillary pressure as a
%                                      function of upscaled saturation,
%                                      current pressure and max. upscaled
%                                      saturation.
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
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('rock2D', [], 'resSat', [], 'pcWG', [], 'p', [], 'rhoW', [], ...
                 'rhoG', []);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.pcWG)
        % We assume a sharp-interface model
        assert(~isempty(opt.resSat));
        sw    = opt.resSat(1);
        sn    = opt.resSat(2);
        h     = Gt.cells.H .* ((S .* (1-sw) - (S_max .* sn)) ./  ...
                                   ((1-sw) .* (1 - sw - sn)));
        h_max = Gt.cells.H .* (S_max ./(1-sw));
        
    else
        % Capillary pressure function provided - we assume a general model
        assert(~isempty(opt.p));
        assert(~isempty(opt.pcWG));
        assert(~isempty(opt.rhoW));
        assert(~isempty(opt.rhoG));
        pc    = opt.pcWG(S, opt.p, 'sGmax', S_max);
        pcmax = opt.pcWG(S_max, opt.p, 'sGmax', S_max);
        drho  = norm(gravity) * (opt.rhoW(opt.p) - opt.rhoG(opt.p));
        h     = pc ./ drho;
        h_max = pcmax ./ drho;
    end
end
