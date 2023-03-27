function debugFz(CO2)

    gravity on;
    g = norm(gravity);
    EOS.rho    = @CO2.rho;
    
    EOS.beta   = @CO2.beta;     
    EOS.beta2  = @(p,t) CO2.rhoDPP(p,t) ./ CO2.rho(p,t);

    EOS.gamma  = @CO2.gamma;        
    EOS.gamma2 = @(p,t) CO2.rhoDTT(p,t) ./ CO2.rho(p,t);        

    EOS.chi    = @(p,t) CO2.rhoDPT(p,t) ./ CO2.rho(p,t);    
    EOS.compressible = 'full';
    
    pI = [1.172834381992800;    
          1.179888067366279;    
          1.187035305199965] * 1e7;

    TI = [3.216031324299658;
          3.218574107533750;
          3.221239130084121] * 1e2;

    h = [0.989584984436863;
         1.046091278527779;
         1.105314001869379] * 1e2; 
    % pI = [1.17; 1.18; 1.17]*1e7; % @@ Artificial case

    temp_grad = 45;

    const_h = false;
    if (const_h)
        h = [1;1;1] * 1e2;
        TI = [1;1;1] * 3.218e2;      % @@ temp. constant on interface
    end
    
    const_temp = false;
    if (const_temp)
        temp_grad = 0;
        TI = [1;1;1] * 3.218e2;      % @@ temp. constant on interface
    end

    rhoI = CO2.rho(pI, TI);    

    [Ieta,  IFpEta,  IFzEta,  Ieta2, Eta, Fp,Fz] = ...
        oldEtaIntegrals(EOS, pI, TI, temp_grad, g);
    
    [NIeta, NIFpEta, NIeta2, NEta, NFp] = etaIntegrals(EOS, pI, TI, temp_grad, g);
    
    dpI = diff(pI);
    dInt = diff(h);
    rhoAvg = avgOf(rhoI);
    
    pT = pI - g * rhoI .* h .* Ieta(-h);
    
    dpT = diff(pT);
    
    dpT2 = avgOf(Fp(-h)) .* dpI + avgOf(Fz{1}(-h)) .* dInt;

    dpT3 = avgOf(NFp(-h)) .* dpI + avgOf( - g * rhoI .*NFp(-h)) .* dInt;
    
    TT = TI - temp_grad * h / 1000; % just for comparison
    rhoT = CO2.rho(pT, TT);         % just for comparison

    keyboard;
end

function avg = avgOf(vec)
    avg = 0.5 * (vec(2:end) + vec(1:end-1));
end



% | Ix |              h e2 |      G |           T_I 1e2 |            P_I e7 |
% |----+-------------------+--------+-------------------+-------------------|
% | 48 | 0.989584984436863 | 317.15 | 3.216031324299658 | 1.172834381992800 |
% | 49 | 1.046091278527779 |        | 3.218574107533750 | 1.179888067366279 |
% | 50 | 1.105314001869379 |        | 3.221239130084121 | 1.187035305199965 |
