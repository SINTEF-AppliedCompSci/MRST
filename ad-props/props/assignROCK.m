function f = assignROCK(f, rock, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    pvtnum = 1;
else
    pvtnum = reg.PVTNUM;
end
cR   = rock(pvtnum, 2);
pRef = rock(pvtnum, 1);

f.cR = cR;
f.pvMultR = @(p)(1 + cR.*(p-pRef));
end

