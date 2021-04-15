function report = wellCalculateProductionComp(state, W, fluid, runspec, time)
%Convert MRST well solution data to ECLIPSE keyword representation.
%
% SYNOPSIS:
%   report = wellCalculateProduction(rSol, wSol, W, fluid, time)
%
% PARAMETERS:
%   rSol  - Reservoir solution structure as defined by a flow and/or a
%           transport solver routine
%
%   wSol  - Well solution data structure.
%
%   W     - Well data structure as defined by function 'addWell'.
%           Report data will be generated for each well in 'W'.
%
%   fluid - Fluid data structure.
%
%   time  - Report time.
%
% RETURNS:
%   report - ECLIPSE keyword representation of the well solution data.
%            A structure array containing the following fields:
%              - TIME -- Report time (== time).
%              - WBHP -- Bottom-hole pressure in each well at TIME.
%              - WVPT -- Total rate (sum of all perforation fluxes).
%              - WOPR -- Total oil production.
%              - WWPR -- Total water production.
%              - WWCT -- Well water-cut.
%
% NOTE:
%   Multiple 'report' structures--e.g., from several time steps--may be
%   concatenated using the 'addToTimeStruct' function.
%
% SEE ALSO:
%   `incompTPFA`, `explicitTransport`, `addWell`, `addToTimeStruct`.

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


   if isempty(W),
      report = struct('TIME', [], 'WBHP', [], 'WVPT', [], ...
                      'WOPR', [], 'WWPR', [], 'WWCT', []);
   else
      phases = {'WATER','OIL','GAS'};
      pl={'W','O','G'};
      for kk=1:numel(phases)
         if(~isfield(runspec,phases{kk}))
            runspec.(phases{kk})=0;
         end
      end
      phcode=size(1,numel(phases));
      for kk=1:numel(phases)
         phcode(kk) = runspec.(phases{kk});
      end
      arg = { 'UniformOutput', false };
      wSol = state.wellSol;
      %mu = fluid.properties(state);
      [c,rho,mu,u] = fluid.pvt(state.pressure,state.z);
      s  = bsxfun(@rdivide,u,sum(u,2));
      %s  = fluid.saturation(state);
      kr = fluid.relperm(s,state);
      %kr = fluid.relperm(s, state);

      mob = bsxfun(@rdivide, kr, mu);
      if(~(size(mob,2)==sum(phcode)))
         error('fluid and phases number do not match');
      end
      if(false) % to compact for hmn to rewrite
         f  = bsxfun(@rdivide, mob, sum(mob, 2));
         fw  = cellfun(@(c) f(c,1), { W.cells }, arg{:});
         rate0 = @(flx, f) sum([bsxfun(@times,flx,f), flx], 1);
         rate1 = @(f,w) rate0(reshape(f, [], 1), reshape(w, [], 1));
         rates = cellfun(rate1, { wSol.flux }, fw, 'UniformOutput', false);
         rates = vertcat(rates{:});  % All 'rates' elements are 1-by-3 DOUBLE.
      else
         f  = bsxfun(@rdivide, mob, sum(mob, 2));
         rates=zeros(numel(W),size(f,2)+1);
         for i=1:numel(W)
            sflux = bsxfun(@times,wSol(i).flux,f(W(i).cells,:));
            % gravity effects not considered
            %press = repmat(wSol.pressure,numel(W(i).cells),1);
            % do upwind weighting assuming only in flux
            press = state.pressure(W(i).cells);
            zz    = state.z(W(i).cells,:);
            %% in effective use, but simple
            [c, rho, mu, u, R,B] = fluid.pvt(press,zz);
            mflux = B\reshape(sflux',[],1);
            mflux = reshape(mflux,size(sflux,2),[])';
            mflux = sum(mflux,1);
            rates(i,:) = [mflux,sum(mflux,2)];
         end
      end


      WBHP = [ wSol.pressure ] .';
      WVPR = rates(:,end);

      report = struct('TIME', time, 'WBHP', WBHP, 'WVPR', WVPR);
      count=0;
      for kk=1:numel(phases)
         if(runspec.(phases{kk})==1)
            count=count+1;
            ff = ['W',pl{kk},'PR'];
            report.(ff)=rates(:,count);
         end
      end
   end
   %eport
end
