function wellSol = getWellSol(W, x, fluid)

mu = fluid.properties();
wellSol(numel(W))=struct;
for i=1:numel(W)
    out = min(x.wellSol(i).flux,0); iout = out<0;
    in  = max(x.wellSol(i).flux,0); iin  = in>0;
    
    swc  = x.s(W(i).cells,:);
    lamc = bsxfun(@rdivide,fluid.relperm(swc), mu); 
    fc   = lamc(:,1)./sum(lamc,2);
    
    sww  = W(i).compi;
    lamw = bsxfun(@rdivide,fluid.relperm(sww), mu);
    fw   = lamw(:,1)./sum(lamw,2);
    
    wellSol(i).name = W(i).name;
    wellSol(i).bhp  = x.wellSol(i).pressure;
    wellSol(i).wcut = iout.*fc + iin.*fw;
    wellSol(i).Sw   = iout.*swc(:,1) + iin.*sww(:,1);
    wellSol(i).qW   = out.*fc + in.*fw;
    wellSol(i).qWs  = sum(wellSol(i).qW);
    wellSol(i).qO   = out.*(1-fc) + in.*(1-fw);
    wellSol(i).qOs  = sum(wellSol(i).qO);
end