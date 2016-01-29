function f = assignSURFROCK(f, surfrock, reg)
ntsfun = numel(reg.SURFINX);
if ntsfun == 1
    surfnum = 1;
else
    surfnum = reg.SURFNUM;
end
f.adsInxSft= surfrock(surfnum, 1);
f.rhoRSft   = surfrock(surfnum, 2);
end
