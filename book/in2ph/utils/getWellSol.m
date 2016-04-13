function wellSol = getWellSol(W, x, fluid)

mu = fluid.properties();
wellSol(numel(W))=struct;
for i=1:numel(W)
    out = min(x.wellSol(i).flux,0); iout = out<0;
    in  = max(x.wellSol(i).flux,0); iin  = in>0;
    lamc = fluid.relperm(x.s(W(i).cells,:))./mu; 
    fc   = lamc(:,1)./sum(lamc,2);
    lamw = fluid.relperm(W(i).compi)./mu;
    fw   = lamw(:,1)./sum(lamw,2);
    wellSol(i).bhp  = x.wellSol(i).pressure;
    wellSol(i).wcut = iout.*fc + iin.*fw;
    wellSol(i).Sw   = iout.*x.s(W(i).cells,1) + iin.*W(i).compi(1);
    wellSol(i).qWs  = sum(out.*fc) + sum(in.*fw);
    wellSol(i).qOs  = sum(out.*(1-fc)) + sum(in.*(1-fw));
end