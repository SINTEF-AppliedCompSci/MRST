function [flag_v, flag_g] = hybridUpwind(Gi, vT)
    nPh = numel(Gi);
    assert(nPh == 2);
    
    flag_v = repmat(vT > 0, 1, nPh);
    dp = Gi{2} - Gi{1};
    d_diff = dp <= 0;
    flag_g = [d_diff, ~d_diff];
end