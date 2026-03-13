function model = Create_kr_history_match(x, model)  
%
% DESCRIPTION: changes the relative permeability in the model struct during
%              the history matching process
%
% SYNOPSIS:
%   model = Create_kr_history_match(x, model)  
%
% PARAMETERS:
%   - x: parameter set for history matching iterations
%   - model: struct on which the history matching process is running
%
% RETURNS:
%   model - struct with the relative permeability changed
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
x = x(:);
row_names = model.history_match.RowNames;
if model.history_match.kr.status 
    if strcmp(model.history_match.kr.type,'MODIFIED-COREY')
        p_s = modified_corey_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = modified_corey_kr(p_s.swc_kr, p_s.sor_kr, p_s.krwSor,...
            p_s.kroSwc, p_s.nW, p_s.nNW);
    elseif strcmp(model.history_match.kr.type,'BROOKS-COREY')
        p_s = brooks_corey_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = brooks_corey_kr(p_s.swc_kr, p_s.sor_kr, p_s.krwSor,...
            p_s.kroSwc, p_s.lambda);
    elseif strcmp(model.history_match.kr.type,'MODIFIED-COREY-MASALMEH')
        p_s = modified_corey_masalmeh_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = modified_corey_masalmeh_kr(p_s.swc_kr, p_s.sor_kr,...
            p_s.krwSor, p_s.kroSwc, p_s.nW, p_s.nNW, p_s.cW, p_s.cNW);           
    elseif strcmp(model.history_match.kr.type,'POINT-BY-POINT')
        sw_kr = sort(model.history_match.kr.Sw_hm);
        sw_kr = sw_kr(:);
        krw = [0; sort(x(1:length(sw_kr)-1))];
        kro = [sort(x(length(sw_kr):2*(length(sw_kr)-1)),'descend'); 0];
        [sw_kr, krw, kro] = kr_correction(sw_kr, krw, kro);
    elseif strcmp(model.history_match.kr.type,'LET')
        p_s = LET_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = LET_kr(p_s.swc_kr, p_s.sor_kr, p_s.krwSor,...
            p_s.kroSwc, p_s.Lw, p_s.Ew, p_s.Tw, p_s.Lnw, p_s.Enw, p_s.Tnw);  
    elseif strcmp(model.history_match.kr.type,'BURDINE')
        sw_pc = model.satfun.sw_pc;
        pc_array = model.satfun.pc;
        p_s = burdine_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = burdine_kr(sw_pc, pc_array, p_s.krwSor, p_s.kroSwc);
    elseif strcmp(model.history_match.kr.type,'VAN-GENUCHTEN-BURDINE')
        p_s = van_genuchten_burdine_kr_handle(model, x, row_names);
        [sw_kr, krw, kro] = van_genuchten_burdine_kr(p_s.swc_kr, p_s.sor_kr,...
             p_s.krwSor, p_s.kroSwc, p_s.n_kr); 
    end
    model.satfun.sw_kr = sw_kr;
    model.satfun.krw = krw;
    model.satfun.kro = kro;

%     figure
%     plot(sw_kr, [krw, kro], 'ks-')
end