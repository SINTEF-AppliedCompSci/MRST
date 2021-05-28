function eqs = addTemperature(eqs,s,dt, state0, T, G,  W, fluid, dp, mob, sF,  bF, bF0, bFqF, pvMult0, pvMult)
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

eqn_num=numel(eqs);
phases=numel(sF);
assert(numel(bFqF)==phases)
assert(numel(mob)==phases)
assert(numel(fluid)==phases+1)
%%
eR = fluid{end}(T);
WT=[W.T]';
eF=cell(phases,1);
eFvF=cell(phases,1);
eFqF=cell(phases,1);
wc=[W.cells]';
for i=1:phases
  upc = (double(dp{i})>=0);
  eF{i}=fluid{i}(T);
  eFvF{i} = s.faceUpstr(upc, eF{i}.*bF{i}.*mob{i}).*s.T.*dp{i};
  %eFw{i}=eF{i}(wc);
  %%
  %eFmobFw(iInx{i})=fluid.eF{i}( WT(iInx{i}) ).*eFmobFw{i};
  %eFmobFw=fluid{i}( WT ).*bFmobFw{i};
  eFqF{i}  = eF{i}(wc).*bFqF{i};
  ind=bFqF{i}<0;
  eFqF{i}(ind)  = fluid{i}( WT(ind) ).*bFqF{i}(ind);
  %eFqF{i}  = fluid{i}( WT ).*bFqF{i};
  %eFqF  = -eFmobFw{i}.*Tw.*(pBHP(perf2well) - pw);% + 0.0*pcOGw + 0.0*g*dzw.*rhoG(wc));
end
dT = s.grad(T);
eRvR  = s.TH.*dT;
eqsn =zeros(G.cells.num,1);
maxE=0;
for i=1:phases
  maxE=max(max(double(fluid{i}(state0.T))),maxE);
  eqsn = eqsn + (s.pv/dt).*( pvMult.*bF{i}.*(eF{i}.*sF{i})- pvMult0.* bF0{i}.*(fluid{i}(state0.T).*state0.s(:,i))) + s.div(eFvF{i});
  eqsn(wc) = eqsn(wc) +  eFqF{i};
end
eqs{eqn_num+1}=eqsn+(s.pv/dt).*( pvMult.*(eR)- pvMult0.*(fluid{end}(state0.T)))+s.div(eRvR);
eqs{eqn_num+1}=eqs{eqn_num+1}./maxE;
end
