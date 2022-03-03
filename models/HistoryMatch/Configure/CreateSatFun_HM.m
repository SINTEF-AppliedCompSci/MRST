function model = CreateSatFun_HM(x, model)  
    x = x(:);
    if model.history_match.kr.status 
        if strcmp(model.history_match.kr.type,'MODIFIED-COREY')
            % Swc ,Sor ,Kwor ,Kowc ,nw ,no 
            Swc = x(1); Sor = x(2); krwSor = x(3); kroSwc = x(4);
            nW = x(5); nNW = x(6);
            [Sw_kr, krw, kro] = modified_corey_kr(Swc, Sor, ...
                krwSor, kroSwc, nW, nNW);
            if model.history_match.pc.status
                if strcmp(model.history_match.pc.type,'BROOKS-COREY')
                    Swc = x(1); Sor = x(2);
                    pd = x(7) * 1e5 ; lambda = x(8);
                    [Sw_pc, pc_array] = brooks_corey_pc(Swc, Sor, pd, lambda);
                elseif strcmp(model.history_match.pc.type,'BROOKS-COREY-SEPARATE-SWC')
                    Swc = x(7); Sor = x(8);
                    pd = x(9) * 1e5 ; lambda = x(10);
                    [Sw_pc, pc_array] = brooks_corey_pc(Swc, Sor, pd, lambda);
                elseif strcmp(model.history_match.pc.type,'VAN-GENUCHTEN')
                    Swc = x(1); Sor = x(2);
                    alpha = x(7) * 1e-5 ; n = x(8);
                    [Sw_pc, pc_array] = van_genuchten(Swc, Sor, alpha, n);
                elseif strcmp(model.history_match.pc.type,'VAN-GENUCHTEN-SEPARATE-SWC')
                    Swc = x(7); Sor = x(8);
                    alpha = x(9) * 1e-5 ; n = x(10);
                    [Sw_pc, pc_array] = van_genuchten(Swc, Sor, alpha, n);
                elseif strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                    %building sw for pc
                    Sw_pc = sort(model.history_match.pc.Sw_hm);
                    Sw_pc = Sw_pc(:);
                    pc_array = sort((x(7:7+length(Sw_pc)-1) .* 1e5),'descend');
                elseif strcmp(model.history_match.pc.type,'SKJAEVELAND')
                    Swc = x(1); Sor = x(2);
                    cwi = x(7) * 1e5; coi = x(8) * 1e5 ; awi = x(9); aoi = x(10);
                    [Sw_pc, pc_array] = skjaeveland_pc(Swc, Sor,...
                                    cwi, coi, awi, aoi);
                elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
                    Swc = x(1); Sor = x(2);
                    cwi = x(7) * 1e5; coi = x(8) * 1e5; ri = x(9) * 1e5; bi = x(10) * 1e5; 
                    Swd = x(11); Sod = x(12);
                    [Sw_pc, pc_array] = modified_skjaeveland_pc...
                        (Swc, Sor, cwi, coi, ri, bi, Swd, Sod);
                end
            end
        elseif strcmp(model.history_match.kr.type,'MODIFIED-COREY-MASALMEH')
            % Swc ,Sor ,Kwor ,Kowc ,nw ,no 
            Swc = x(1); Sor = x(2); 
            krwSor = x(3); kroSwc = x(4);
            nW = x(5); nNW = x(6); cW = x(7); cNW = x(8);
            [Sw_kr, krw, kro] = modified_corey_masalmeh_kr(Swc, Sor, ...
                krwSor, kroSwc, nW, nNW, cW, cNW);
            if model.history_match.pc.status
                if strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND-MASALMEH')
                    Swc = x(1); Sor = x(2);
                    cwi = x(9) * 1e5; coi = x(10) * 1e5 ; awi = x(11); aoi = x(12);
                    sw_cutoff = x(13); b = x(14) * 1e5;
                    [Sw_pc, pc_array] = modified_skjaeveland_masalmeh_pc...
                        (Swc, Sor, cwi, coi, awi, aoi, sw_cutoff, b);
                elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND-MASALMEH-SEPARATE-SWC')
                    Swc = x(9); Sor = x(10);
                    cwi = x(11) * 1e5; coi = x(12) * 1e5 ; awi = x(13); aoi = x(14);
                    sw_cutoff = x(15); b = x(16) * 1e5;
                    [Sw_pc, pc_array] = modified_skjaeveland_masalmeh_pc...
                        (Swc, Sor, cwi, coi, awi, aoi, sw_cutoff, b);
                end
            end            
        elseif strcmp(model.history_match.kr.type,'POINT-BY-POINT')
            Sw_kr = sort(model.history_match.kr.Sw_hm);
            Sw_kr = Sw_kr(:);
            krw = [0; sort(x(1:length(Sw_kr)-1))];
            kro = [sort(x(length(Sw_kr):2*(length(Sw_kr)-1)),'descend'); 0];
            [Sw_kr, krw, kro] = kr_correction(Sw_kr, krw, kro);
            if model.history_match.pc.status
                if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                    Sw_pc = sort(model.history_match.pc.Sw_hm);
                    Sw_pc = Sw_pc(:);
                    Sw_kr_hm = sort(model.history_match.kr.Sw_hm);
                    Sw_kr_hm = Sw_kr_hm(:); len_sw_kr = numel(Sw_kr_hm);
                    pc_array = sort((x(len_sw_kr*2-1:len_sw_kr*2+length(Sw_pc)-2) .* 1e5),'descend');
                end
            end
        elseif strcmp(model.history_match.kr.type,'LET')
            Swc = x(1); Sor = x(2); krwSor = x(3); kroSwc = x(4);
            Lw = x(5); Lnw = x(6); Ew = x(7); Enw = x(8); Tw = x(9); Tnw = x(10);
            [Sw_kr, krw, kro] = LET_kr(Swc, Sor, krwSor, kroSwc,...
                Lw, Ew, Tw, Lnw, Enw, Tnw);  
            if model.history_match.pc.status
                if strcmp(model.history_match.pc.type,'LET-DRAINAGE')
                    Swc = x(1); max_pc = x(11) * 1e5;
                    entry_pc = x(12) * 1e5;
                    entry_multiplier = x(13);
                    forced_multiplier = x(14);
                    L_entry = x(15); 
                    E_entry = x(16); 
                    T_entry = x(17);
                    L_forced = x(18);
                    E_forced = x(19);
                    T_forced = x(20);
                    [Sw_pc, pc_array] = LET_drainage_pc(Swc, entry_multiplier,...
                        forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
                        E_forced, T_forced);
                elseif strcmp(model.history_match.pc.type,'LET-DRAINAGE-SEPARATE-SWC')
                    Swc = x(11); max_pc = x(12) * 1e5;
                    entry_pc = x(13) * 1e5;
                    entry_multiplier = x(14);
                    forced_multiplier = x(15);
                    L_entry = x(16); 
                    E_entry = x(17); 
                    T_entry = x(18);
                    L_forced = x(19);
                    E_forced = x(20);
                    T_forced = x(21);
                    [Sw_pc, pc_array] = LET_drainage_pc(Swc, entry_multiplier,...
                        forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
                        E_forced, T_forced);
                elseif strcmp(model.history_match.pc.type,'LET-IMBIBITION')
                    Swc = x(1); Sor = x(2); 
                    sw_pc0 = x(11);
                    max_pc = x(12) * 1e5;
                    min_pc = x(13) * 1e5;
                    spontaneous_multiplier = x(14);
                    forced_multiplier = x(15);
                    L_spont = x(16); E_spont = x(17); T_spont = x(18);
                    L_forced = x(19); E_forced = x(20); T_forced = x(21);
                    [Sw_pc, pc_array] = LET_imbibition_pc(Swc, Sor, spontaneous_multiplier,...
                        forced_multiplier, min_pc, max_pc, L_spont, E_spont, T_spont, L_forced, ...
                        E_forced, T_forced, sw_pc0);
                elseif strcmp(model.history_match.pc.type,'BROOKS-COREY')
                    Swc = x(1); Sor = x(2);
                    pd = x(11) * 1e5 ; lambda = x(12);
                    [Sw_pc, pc_array] = brooks_corey_pc(Swc, Sor, pd, lambda);
                elseif strcmp(model.history_match.pc.type,'POINT-BY-POINT')
                    %building sw for pc
                    Sw_pc = sort(model.history_match.pc.Sw_hm);
                    Sw_pc = Sw_pc(:);
                    pc_array = sort((x(11:11+length(Sw_pc)-1) .* 1e5),'descend');
                elseif strcmp(model.history_match.pc.type,'SKJAEVELAND')
                    Swc = x(1); Sor = x(2);
                    cwi = x(11) * 1e5; coi = x(12) * 1e5; awi = x(13); aoi = x(14);
                    [Sw_pc, pc_array] = skjaeveland_pc(Swc, Sor,...
                        cwi, coi, awi, aoi);
                elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
                    Swc = x(1); Sor = x(2);
                    cwi = x(11) * 1e5; coi = x(12) * 1e5; ri = x(13) * 1e5;
                    bi = x(14) * 1e5; Swd = x(15); Sod = x(16);
                    [Sw_pc, pc_array] = modified_skjaeveland_pc...
                        (Swc, Sor, cwi, coi, ri, bi, Swd, Sod);
                end
            end
        end
    end
    if model.history_match.pc.status && not(model.history_match.kr.status)
        if strcmp(model.history_match.pc.type,'POINT-BY-POINT')
            %building sw for pc
            Sw_pc = sort(model.history_match.pc.Sw_hm);
            Sw_pc = Sw_pc(:);
            pc_array = sort((x(1:length(Sw_pc)) .* 1e5),'descend');
        elseif strcmp(model.history_match.pc.type,'BROOKS-COREY')
            Swc = x(1); Sor = x(2);
            pd = x(3) * 1e5 ; lambda = x(4);
            [Sw_pc, pc_array] = brooks_corey_pc(Swc, Sor, pd, lambda);
        elseif strcmp(model.history_match.pc.type,'SKJAEVELAND')
            Swc = x(1); Sor = x(2);
            cwi = x(3) * 1e5; coi = x(4) * 1e5; awi = x(5);  aoi = x(6);
            [Sw_pc, pc_array] = skjaeveland_pc(Swc, Sor,...
                cwi, coi, awi, aoi);
        elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
            Swc = x(1); Sor = x(2);
            cwi = x(3) * 1e5; coi = x(4) * 1e5 ; ri = x(5) * 1e5; bi = x(6) * 1e5;
            Swd = x(7); Sod = x(8);
            [Sw_pc, pc_array] = modified_skjaeveland_pc...
                (Swc, Sor, cwi, coi, ri, bi, Swd, Sod);
        elseif strcmp(model.history_match.pc.type,'LET-IMBIBITION')
            Swc = x(1); Sor = x(2); 
            sw_pc0 = x(3);
            max_pc = x(4) * 1e5;
            min_pc = x(5) * 1e5;
            spontaneous_multiplier = x(6);
            forced_multiplier = x(7);
            L_spont = x(8); E_spont = x(9); T_spont = x(10);
            L_forced = x(11); E_forced = x(12); T_forced = x(13);
            [Sw_pc, pc_array] = LET_imbibition_pc(Swc, Sor, spontaneous_multiplier,...
                forced_multiplier, min_pc, max_pc, L_spont, E_spont, T_spont, L_forced, ...
                E_forced, T_forced, sw_pc0);
        elseif strcmp(model.history_match.pc.type,'LET-DRAINAGE')
            Swc = x(1); max_pc = x(2) * 1e5;
            entry_pc = x(3) * 1e5;
            entry_multiplier = x(4);
            forced_multiplier = x(5);
            L_entry = x(6); 
            E_entry = x(7); 
            T_entry = x(8);
            L_forced = x(9);
            E_forced = x(10);
            T_forced = x(10);
            [Sw_pc, pc_array] = LET_drainage_pc(Swc, entry_multiplier,...
                forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
                E_forced, T_forced);
        end
    end

    if model.history_match.kr.status
        model.satfun.sw_kr = Sw_kr;
        model.satfun.krw = krw;
        model.satfun.kro = kro;
    end
    if model.history_match.pc.status
        model.satfun.sw_pc = Sw_pc;
        model.satfun.pc = pc_array;
    end
    
end