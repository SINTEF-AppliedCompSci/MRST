function eqs = addTemperature(eqs,s,dt, state0, T, G,  W, fluid, dp, mob, sF,  bF, bF0, bFqF, pvMult0, pvMult)
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
