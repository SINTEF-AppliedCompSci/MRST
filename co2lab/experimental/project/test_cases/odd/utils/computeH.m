function h = computeH(EOS, P_top, T_top, Tgrad, rhoW, g_cos_t, M)

    h0 = M ./ EOS(P_top, T_top);  % result should be a scalar
    
    % since 'fzero' must be called with scalar arguments, we must make a loop
    % here.
    
    for i = 1:numel(M)
        if (M(i) > 0)
            h(i) = fzero(@(x)residual(x,i), h0(i));
        else
            h(i) = 0;
        end
    end
        
    function res = residual(height, i)
        P_int   = P_top(i) + rhoW * g_cos_t * height;
        T_int   = T_top(i) + Tgrad * height / 1000;
        rho_int = EOS(P_int, T_int);
        Ieta    = etaIntegrals(EOS, P_int, T_int, Tgrad, g_cos_t);
        res     = height * Ieta(-height) * rho_int - M(i);
    end

end