function model = CreateKr(model)
%
% DESCRIPTION: adds relative permeability table to the model
%
% SYNOPSIS:
%   model = CreateKr(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment: saturation functions used for forward modeling
%
% RETURNS:
%   model - adds the the following fields to the struct:
%   - satfun: struct with relative permeability information inside
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
kr_struct = model.experiment.satfun.kr;
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
            kr_struct.krwSor, kr_struct.kroSwc, kr_struct.nW, kr_struct.nNW,...
            kr_struct.cW, kr_struct.cNW);
    case 'BROOKS-COREY'
        lambda = kr_struct.lambda;
        [Sw, krw, kro] = brooks_corey_kr(kr_struct.Swc, kr_struct.Sor, ...
            kr_struct.krwSor, kr_struct.kroSwc, lambda);
    case 'BURDINE'
        Sw_pc = model.satfun.sw_pc;
        pc_array = model.satfun.pc;
        [Sw, krw, kro] = burdine_kr(Sw_pc, pc_array, kr_struct.krwSor,...
            kr_struct.kroSwc);
    case 'VAN-GENUCHTEN-BURDINE'
        n = kr_struct.n;
        [Sw, krw, kro] = van_genuchten_burdine_kr(kr_struct.Swc, ....
            kr_struct.Sor, kr_struct.krwSor, kr_struct.kroSwc, n);
    case 'LET'
        [Sw, krw, kro] = LET_kr(kr_struct.Swc, kr_struct.Sor, ...
            kr_struct.krwSor, kr_struct.kroSwc, kr_struct.Lw, ...
            kr_struct.Ew, kr_struct.Tw, kr_struct.Lnw, kr_struct.Enw,...
            kr_struct.Tnw);  
end

model.satfun.sw_kr = Sw;
model.satfun.krw = krw;
model.satfun.kro = kro;
