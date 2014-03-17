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
%                - pcOG(S, p, S_max) - upscaled capillary pressure as a
%                                      funciton of upscaled saturation,
%                                      current pressure and max. upscaled
%                                      saturation.
%                - rhoO(p)           - density of water [oil] phase, as a
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
%
% EXAMPLE:
%
% SEE ALSO:
%
    opt = struct('rock2D', [], 'resSat', [], 'pcOG', [], 'p', [], 'rhoO', [], ...
                 'rhoG', []);
    opt = merge_options(opt, varargin{:});
    
    if isempty(opt.pcOG)
        % We assume a sharp-interface model
        assert(~isempty(opt.resSat));
        sw    = opt.resSat(1);
        sn    = opt.resSat(2);
        h     = Gt.cells.H .* (S .* (1-sw) - (S_max .* sn)) ./  ...
                                   ((1-sw) .* (1 - sw - sn));
        h_max = Gt.cells.H .* S_max ./(1-sw);
        
    else
        % Capillary pressure function provided - we assume a general model
        assert(~isempty(opt.p));
        assert(~isempty(p));
        assert(~isempty(rhoO));
        assert(~isempty(rhoG));
        pc    = opt.pcOG(S, p, 'sGmax', S_max);
        pcmax = opt.pcOG(S_max, p, 'sGmax', S_max);
        drho  = norm(gravity) * (rhoO(p) - rhoG(p));
        h     = pc ./ drho;
        h_max = pcmax ./ drho;
    end
end
