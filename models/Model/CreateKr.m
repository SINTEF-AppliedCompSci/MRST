function model = CreateKr(model)
% <keywords>
%
% Purpose : create a relative permeability table
%
% Syntax :
%   model = CreateKr(model)
%
% Input Parameters :
%   model: struct output from createrock function
%
% Return Parameters :
%   model: struct containing the realtive permeability table
%
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
    kr_struct = model.experiment.satfun.kr;
    pc_struct = model.experiment.satfun.pc;
    switch string(kr_struct.type)
        case 'TABLE'
            [Sw, ia]  = unique(kr_struct.table.(1)); 
            kr_struct.table = kr_struct.table(ia,:);
            krw = kr_struct.table.(2); 
            kro = kr_struct.table.(3);    
            kr_struct.Swc    = min(Sw); 
            kr_struct.Sor    = 1 - max(Sw);
            kr_struct.krwSor = max(krw);
            kr_struct.kroSwc = max(kro); 
            model.experiment.satfun.kr = kr_struct;
            [Sw, krw, kro] = kr_correction(Sw, krw, kro);
        case 'MODIFIED-COREY'
            [Sw, krw, kro] = modified_corey_kr(kr_struct.Swc, kr_struct.Sor, ...
                kr_struct.krwSor, kr_struct.kroSwc, kr_struct.nW, kr_struct.nNW);
        case 'MODIFIED-COREY-MASALMEH'
            [Sw, krw, kro] = modified_corey_masalmeh_kr(kr_struct.Swc, kr_struct.Sor, ...
                kr_struct.krwSor, kr_struct.kroSwc, kr_struct.nW, kr_struct.nNW, kr_struct.cW, kr_struct.cNW);
        case 'BROOKS-COREY'
            lambda = kr_struct.lambda;
            [Sw, krw, kro] = brooks_corey_kr(kr_struct.Swc, kr_struct.Sor, kr_struct.krwSor,...
                kr_struct.kroSwc, lambda);
        case 'BURDINE'
            Sw_pc = model.satfun.sw_pc;
            pc_array = model.satfun.pc;
            [Sw, krw, kro] = burdine_kr(Sw_pc, pc_array, kr_struct.krwSor, kr_struct.kroSwc);
        case 'VAN-GENUCHTEN-BURDINE'
            n = kr_struct.n;
            [Sw, krw, kro] = van_genuchten_burdine_kr(kr_struct.Swc, ....
                kr_struct.Sor, kr_struct.krwSor, kr_struct.kroSwc, n);
        case 'LET'
            [Sw, krw, kro] = LET_kr(kr_struct.Swc, kr_struct.Sor, kr_struct.krwSor, kr_struct.kroSwc,...
                kr_struct.Lw, kr_struct.Ew, kr_struct.Tw, kr_struct.Lnw, kr_struct.Enw, kr_struct.Tnw);  
    end
    
    model.satfun.sw_kr = Sw;
    model.satfun.krw = krw;
    model.satfun.kro = kro;
end
