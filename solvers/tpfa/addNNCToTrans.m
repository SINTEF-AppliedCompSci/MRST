function T = addNNCToTrans(T, trans_nnc)
% Transmissibility to half-transmissibility
    T = [T; 2*trans_nnc; 2*trans_nnc];
end