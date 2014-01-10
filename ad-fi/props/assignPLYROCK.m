function f = assignPLYROCK(f, plyrock, reg)
ntsfun = numel(reg.SATINX);
if ntsfun == 1
    satnum = 1;
else
    satnum = reg.SATNUM;
end
f.dps    = plyrock(satnum, 1);
f.rrf    = plyrock(satnum, 2);
f.rhoR   = plyrock(satnum, 3);
f.adsInx = plyrock(satnum, 4);
f.adsMax   = plyrock(satnum, 5);
end
