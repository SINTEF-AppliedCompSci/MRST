function f = assignSURFROCK(f, surfrock, reg)
ntsfun = numel(reg.SATINX);
if ntsfun == 1
    satnum = 1;
else
    satnum = reg.SATNUM;
end
f.adsSurfInx = surfrock(satnum, 1);
f.rhoSurfR   = surfrock(satnum, 2);
end
