function [sol, withinLims] = updateControls(W, sol, pBH, q_s, model, varargin)
%Undocumented Utility Function

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

   opt = struct('Verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   if isempty(W)||isempty(W(1).lims)
      withinLims = true(numel(sol),1);
   else

      nwells     = numel(sol);
      withinLims = true(nwells,1);

      pBH   = value(pBH);
      q_s   = cell2mat( cellfun(@value, q_s, 'UniformOutput', false) );
      q_s = getRatesWOG(q_s, model);

      for wnr = 1:numel(sol)
         lims = W(wnr).lims;
         if ~isfield(lims, 'vrat')
            lims.vrat = -Inf;
         end
         pBHw  = pBH(wnr);
         q_sw  = q_s(wnr,:);
         qt_sw = sum(q_sw);
         if ~isnumeric(W(wnr).lims)
            if sol(wnr).sign > 0   % injector
               modes   = {'bhp', 'rate'};
               flags = [pBHw > lims.bhp, qt_sw > lims.rate];
            else            % producer
               modes   = {'bhp', 'orat', 'lrat', 'grat', 'wrat', 'vrat'};
               flags = [pBHw           < lims.bhp,  ...
                        q_sw(2)         < lims.orat, ...
                        q_sw(1)+q_sw(2) < lims.lrat, ...
                        q_sw(3)         < lims.grat, ...
                        q_sw(1)         < lims.wrat, ...
                        qt_sw           < lims.vrat];
            end
         else
            modes = {};
            flags = false(numel(sol), 1);
            assert(isinf(lims))
         end
         %limits we need to check (all others than w.type):
         chkInx = ~strcmp(sol(wnr).type, modes);
         vltInx = find(flags(chkInx), 1);
         if ~isempty(vltInx)
            withinLims(wnr) = false;
            modes  = modes(chkInx);
            switchMode = modes{vltInx};
            if opt.Verbose
               fprintf('Well %s: Control mode changed from %s to %s.\n', sol(wnr).name, sol(wnr).type, ...
                       switchMode);
            end
            sol(wnr).type = switchMode;
            sol(wnr).val  = lims.(switchMode);
         end
      end
   end
end

function q_s = getRatesWOG(q_s, model)
   switch model
     case 'OW'
       q_s = [q_s, zeros(size(q_s,1), 1)];
     case {'3P', 'BO', 'VO'}
       % nothing
     otherwise
       error(['Unknown model: ', model]);
   end
end








