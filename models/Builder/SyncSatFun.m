function [Sw,krw,kro,pc] = SyncSatFun(Sw,krw,kro,Sw_pc,pc)
    % make kr to the same length of pc
    if (length(Sw_pc) > length(Sw))
        krw = interp1(Sw,krw,Sw_pc);
        kro = interp1(Sw,kro,Sw_pc);
        Sw  = Sw_pc;
    end
    % make pc to the same length of kr
    if (length(Sw) > length(Sw_pc))
        pc = interp1(Sw_pc,pc,Sw);
    end
end