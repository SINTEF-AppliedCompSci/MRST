function [krW, krO] = relPermHelper(sW, krW, krO)
    krW = krW(sW);
    krO = krO(1 - sW);
end
