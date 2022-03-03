function model = CreatePc(model)
% <keywords>
%
% Purpose : create a capillary pressure saturation function table
%
% Syntax :
%   model = CreatePc(model)
%
% Input Parameters :
%   model: struct output from the CreateKr function
%
% Return Parameters :
%   model: struct containing the capillary pressure table
% Description :
%
% Author : 
%    Siroos Azizmohammadi
%    Omidreza Amrollahinasab
%
% History :
% \change{1.0}{09-Nov-2021}{Original}
%
% --------------------------------------------------
% (c) 2021, Siroos Azizmohammadi,
% Omidreza Amrollahinasab
% Chair of Reservoir Engineering, University of Leoben, Austria
% email: info@dpe.ac.at
% url: dpe.ac.at
% --------------------------------------------------
%
%%
    pc_struct = model.experiment.satfun.pc;
    kr_struct = model.experiment.satfun.kr;
    switch string(pc_struct.type)
        case 'TABLE'
            [Sw_pc,ia] = unique(pc_struct.table.(1));
            pc_struct.table = pc_struct.table(ia,:);
            pc_array = pc_struct.table.(2);
        case 'BROOKS-COREY'
            lambda = pc_struct.lambda;
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            [Sw_pc, pc_array] = brooks_corey_pc(kr_struct.Swc, kr_struct.Sor, pc_struct.pd, lambda);
        case 'BROOKS-COREY-SEPARATE-SWC'
            lambda = pc_struct.lambda;
            [Sw_pc, pc_array] = brooks_corey_pc(pc_struct.Swc, pc_struct.Sor, pc_struct.pd, lambda);
        case 'VAN-GENUCHTEN'
            n = pc_struct.n;
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            [Sw_pc, pc_array] = van_genuchten(kr_struct.Swc, kr_struct.Sor, pc_struct.alpha, n);
        case 'VAN-GENUCHTEN-SEPARATE-SWC'
            n = pc_struct.n;
            [Sw_pc, pc_array] = van_genuchten(pc_struct.Swc, pc_struct.Sor, pc_struct.alpha, n);
        case 'SKJAEVELAND'            
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            [Sw_pc, pc_array] = skjaeveland_pc(kr_struct.Swc, kr_struct.Sor,...
                pc_struct.cwi, pc_struct.coi, pc_struct.awi, pc_struct.aoi);
        case 'MODIFIED-SKJAEVELAND-MASALMEH'   
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            cwi = pc_struct.cwi; awi = pc_struct.awi;
            coi = pc_struct.coi; aoi = pc_struct.aoi;
            sw_cutoff = pc_struct.sw_cutoff; b = pc_struct.b;
            [Sw_pc, pc_array] = modified_skjaeveland_masalmeh_pc...
                (kr_struct.Swc, kr_struct.Sor, cwi, coi, awi, aoi, sw_cutoff, b);
        case 'MODIFIED-SKJAEVELAND-MASALMEH-SEPARATE-SWC'            
            [Sw_pc, pc_array] = modified_skjaeveland_masalmeh_pc...
                (pc_struct.Swc, pc_struct.Sor, pc_struct.cwi, ...
                pc_struct.coi, pc_struct.awi, pc_struct.aoi,...
                pc_struct.sw_cutoff, pc_struct.b);
        case 'MODIFIED-SKJAEVELAND'
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            [Sw_pc, pc_array] = modified_skjaeveland_pc...
                (kr_struct.Swc, kr_struct.Sor, pc_struct.cwi, pc_struct.coi,...
                pc_struct.ri, pc_struct.bi, pc_struct.Swd, pc_struct.Sod);
        case 'ZERO'
            nPoints = 50;
            Sw_pc  = linspace(0, 1, nPoints)';
            pc_array = zeros(length(Sw_pc),1);     
        case 'LET-IMBIBITION'
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            spontaneous_multiplier = pc_struct.spontaneous_multiplier;
            forced_multiplier = pc_struct.forced_multiplier;
            min_pc = pc_struct.min_pc; max_pc = pc_struct.max_pc;
            L_spont = pc_struct.L_spont; E_spont = pc_struct.E_spont;
            T_spont = pc_struct.T_spont; L_forced = pc_struct.L_forced;
            E_forced = pc_struct.E_forced; T_forced = pc_struct.T_forced;
            sw_pc0 = pc_struct.sw_pc0;
            [Sw_pc, pc_array] = LET_imbibition_pc(kr_struct.Swc, kr_struct.Sor, spontaneous_multiplier,...
                forced_multiplier, min_pc, max_pc, L_spont, E_spont, T_spont, L_forced, ...
                E_forced, T_forced, sw_pc0);
        case 'LET-DRAINAGE'
            fprintf('calculating pc table: using Swc and Sor from kr model...\n')
            max_pc = pc_struct.max_pc;
            entry_pc = pc_struct.entry_pc;
            entry_multiplier = pc_struct.entry_multiplier;
            forced_multiplier = pc_struct.forced_multiplier;
            L_entry = pc_struct.L_entry; E_entry = pc_struct.E_entry; T_entry = pc_struct.T_entry;
            L_forced = pc_struct.L_forced; E_forced = pc_struct.E_forced; T_forced = pc_struct.T_forced;
            Swc = kr_struct.Swc;
            [Sw_pc, pc_array] = LET_drainage_pc(Swc, entry_multiplier,...
                forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
                E_forced, T_forced);
        case 'LET-DRAINAGE-SEPARATE-SWC'
            max_pc = pc_struct.max_pc;
            entry_pc = pc_struct.entry_pc;
            entry_multiplier = pc_struct.entry_multiplier;
            forced_multiplier = pc_struct.forced_multiplier;
            L_entry = pc_struct.L_entry; E_entry = pc_struct.E_entry; T_entry = pc_struct.T_entry;
            L_forced = pc_struct.L_forced; E_forced = pc_struct.E_forced; T_forced = pc_struct.T_forced;
            Swc = pc_struct.Swc;
            [Sw_pc, pc_array] = LET_drainage_pc(Swc, entry_multiplier,...
                forced_multiplier, entry_pc, max_pc, L_entry, E_entry, T_entry, L_forced, ...
                E_forced, T_forced);
    end
    
    model.satfun.sw_pc = Sw_pc;
    model.satfun.pc = pc_array;
end
