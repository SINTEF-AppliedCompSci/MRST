function model = Create_pc_history_match(x, model, params_no_kr)  
%
% DESCRIPTION: changes the capillary pressure in the model struct during
%              the history matching process
%
% SYNOPSIS:
%   model = Create_pc_history_match(x, model, params_no_kr)  
%
% PARAMETERS:
%   - x: parameter set for history matching iterations
%   - model: struct on which the history matching process is running
%   - params_no_kr: number of the parameters for the relative permeability
%
% RETURNS:
%   model - struct with the capillary pressure changed
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
if model.history_match.pc.status
    if strcmp(model.history_match.pc.type,'BROOKS-COREY')
        p_s = brooks_corey_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = brooks_corey_pc(p_s.swc_pc, p_s.sor_pc, p_s.pd, p_s.lambda_pc);
    elseif strcmp(model.history_match.pc.type,'VAN-GENUCHTEN')
        p_s = van_genuchten_handle(model, x, row_names);
        [sw_pc, pc_array] = van_genuchten(p_s.swc_pc, p_s.sor_pc, p_s.alpha, p_s.n_pc);
    elseif strcmp(model.history_match.pc.type,'POINT-BY-POINT')
        if model.history_match.kr.status && strcmp(model.history_match.kr.type,'POINT-BY-POINT')
            sw_pc = sort(model.history_match.pc.Sw_hm);
            sw_pc = sw_pc(:);
            Sw_kr_hm = sort(model.history_match.kr.Sw_hm);
            Sw_kr_hm = Sw_kr_hm(:); len_sw_kr = numel(Sw_kr_hm);
            pc_array = sort((x(len_sw_kr*2-1:len_sw_kr*2+length(sw_pc)-2) .* Convert('bar')),'descend');
        else
            sw_pc = sort(model.history_match.pc.Sw_hm);
            sw_pc = sw_pc(:);
            pc_array = sort((x(params_no_kr+1:params_no_kr+1+length(sw_pc)-1) .* Convert('bar')),'descend');
        end
    elseif strcmp(model.history_match.pc.type,'SKJAEVELAND')
        p_s = skjaeveland_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = skjaeveland_pc(p_s.swc_pc, p_s.sor_pc, p_s.cwi, p_s.coi, p_s.awi, p_s.aoi);
    elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND')
        p_s = modified_skjaeveland_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = modified_skjaeveland_pc(p_s.swc_pc, p_s.sor_pc, p_s.cwi, p_s.coi, p_s.ri, p_s.bi, p_s.swd, p_s.sod);
    elseif strcmp(model.history_match.pc.type,'MODIFIED-SKJAEVELAND-MASALMEH')
        p_s = modified_skjaeveland_masalmeh_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = modified_skjaeveland_masalmeh_pc(...
            p_s.swc_pc, p_s.sor_pc, p_s.cwi, p_s.coi, p_s.awi, p_s.aoi, p_s.sw_cutoff, p_s.b);
    elseif strcmp(model.history_match.pc.type,'LET-DRAINAGE')
        p_s = LET_drainage_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = LET_drainage_pc(p_s.swc_pc, p_s.entry_multiplier,...
            p_s.forced_multiplier, p_s.entry_pc, p_s.max_pc, p_s.L_entry, p_s.E_entry, p_s.T_entry, p_s.L_forced, ...
            p_s.E_forced, p_s.T_forced);
    elseif strcmp(model.history_match.pc.type,'LET-IMBIBITION')
        p_s = LET_imbibition_pc_handle(model, x, row_names);
        [sw_pc, pc_array] = LET_imbibition_pc(p_s.swc_pc, p_s.sor_pc, p_s.spontaneous_multiplier,...
            p_s.forced_multiplier, p_s.min_pc, p_s.max_pc, p_s.L_spont, p_s.E_spont, p_s.T_spont, p_s.L_forced, ...
            p_s.E_forced, p_s.T_forced, p_s.sw_pc0);
    end
    model.satfun.sw_pc = sw_pc;
    model.satfun.pc = pc_array;

%     figure
%     plot(sw_pc, pc_array / Convert('bar'), 'ks-')
end   