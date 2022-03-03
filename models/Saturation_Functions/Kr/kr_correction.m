function [Sw, krw, kro] = kr_correction(Sw, krw, kro)

    % adding tale ends
    if not(Sw(1) == 0)
        Sw = sort(Sw); krw = sort(krw);
        kro = sort(kro,'descend');
        Sw = [0;Sw]; krw = [0;krw];
        kro = [max(kro);kro];
    end
    if not(Sw(end) == 1)
        Sw = sort(Sw); krw = sort(krw);
        kro = sort(kro,'descend');
        Sw = [Sw;1]; krw = [krw;max(krw)];
        kro = [kro;0];
    end
    
    % check for invalid entries
    krw(or(Sw>1,Sw<0)) = [];
    kro(or(Sw>1,Sw<0)) = [];
    Sw(or(Sw>1,Sw<0)) = [];
end
