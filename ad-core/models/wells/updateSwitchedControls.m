function [W, wells_shut] = updateSwitchedControls(wellSol, W, varargin)
%Update active well controls based on well solution structure
%
% SYNOPSIS:
%   [W, wells_shut] = updateSwitchedControls(sol, W)
%
% DESCRIPTION:
%   Updates the well control structure based on the solution structure.
%   Typically, if a well has changed controls during a timestep, we want to
%   modify the controls for the next timestep accordingly. The purpose of
%   this function is twofold:
%           - Update the well controls based on limits
%           - Shut down any wells which have changed from injecting to
%           producing or vice versa.
%
%   Note that either feature can be enabled/disabled via keyword arguments.
%
% REQUIRED PARAMETERS:
%   wellSol - Well sol struct, corresponding to the solution of a problem
%             using the current set of wells. The wellSol will contain
%             information about how the controls have changed during the
%             timestep.
%
%   W       - Wells used to produce wellSol.
%
% OPTIONAL PARAMETERS:
%   allowWellSignChange - Boolean indicating if wells are allowed to change
%                         between injection and production. Disabled by
%                         default (changing wells will be shut down until
%                         re-enabled by new controls).
%
%  allowControlSwitching- Boolean indicating if wells are allowed to change
%                         controls based om limits (W.lims). Default true.
%
% RETURNS:
%   W       - Wells with controls updated based on wellSol.
%
%   wells_shut - Boolean indicating if any wells were shut down due to
%                changing from injection to production or vice versa,
%                dependent on allowWellSignChange keyword argument.
%
% SEE ALSO:
%   WellModel

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

   wells_shut = false;
   if isempty(W) || isempty(wellSol)
      return
   end
   
   if isfield(W, 'status')
       active = vertcat(W.status);
   else 
       active = true(size(W));
   end
   
   if ~any(active)
       return
   end
   
   W_tmp = W(active);
   wellSol   = wellSol(active);
   
   opt = struct('allowWellSignChange', false, ...
                'allowControlSwitching', true, ...
                'Verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   % Check if producers are becoming injectors and vice versa. The indexes
   % of such wells are stored in inx.
   wsg = vertcat(W_tmp(:).sign);
   ssg = sign(getTotalRate(wellSol));
   inx = wsg ~= ssg;

   % A well can be set to zero rate without beeing shut down. We update inx
   % to take into account this fact.
   wval = vertcat(W_tmp(:).val);
   wtype = {W_tmp(:).type}';
   inx = inx & ~strcmp('bhp', wtype) & (wval ~= 0);

   inx = find(inx);
   if ~opt.allowWellSignChange && opt.allowControlSwitching
      if any(inx)
         wells_shut = true;
         ostring = 'Wells ';

         for k = 1:numel(inx)
            W_tmp(inx(k)).status = false;
            ostring = [ostring,  W_tmp(inx(k)).name, ', '];
         end

         if opt.Verbose,
            fprintf([ostring(1:end - 2), ' shut down.\n']);
         end
      end
   end

   % Check if well-controls have been switch, if so, update W
   inx = find(~arrayfun(@(x,y)strcmp(x.type,y.type), W_tmp(:), wellSol(:)));
   for k = 1:numel(inx)
      fromTp = W_tmp(inx(k)).type;
      toTp   = wellSol(inx(k)).type;

      if opt.Verbose
         fprintf(['Well ', W_tmp(inx(k)).name,       ...
                  ' has switched from ', fromTp, ...
                  ' to ', toTp, '.\n']);
      end

      W_tmp(inx(k)).type = toTp;
      W_tmp(inx(k)).val  = wellSol(inx(k)).val;
   end
   W(active) = W_tmp;
end

%--------------------------------------------------------------------------

function qt = getTotalRate(sol)
   ns = numel(sol);
   qt       = zeros([ns, 1]);
   if ns == 0
       return
   end
   typelist = {'qWs', 'qOs', 'qGs'};
   types    = typelist(isfield(sol(1), typelist));
   for w = 1:ns
      for t = reshape(types, 1, []),
         x = sol(w).(t{1});
         if ~isempty(x),
            qt(w) = qt(w) + x;
         end
      end
   end
end
