function f = assignRSCONSTT(f, rsconstt, reg)
ntpvt = numel(reg.PVTINX);
if ntpvt == 1
    f.rs   = rsconstt(1,1);
    f.pBub = rsconstt(1,2);
else
    f.rs   = rsconstt(reg.PVTNUM,1);
    f.pBub = rsconstt(reg.PVTNUM,2);
end
