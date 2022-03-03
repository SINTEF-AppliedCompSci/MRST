function pc_struct = read_pc_parameters_from_settings(pc_struct, line)

if(strcmp(pc_struct.type,strcat('TABLE')))
    pc_struct.fullFile = line{3};
    [filePath,fileName,fileExt] = fileparts(pc_struct.fullFile);
    pc_struct.filePath = filePath;
    pc_struct.fileName = strcat(fileName,fileExt);
    [~,table] = ImportTable(pc_struct.fullFile);                     
%     [~,inputTable,~,table] = ImportTable2(pc.fullFile);
%     pc.inputTable = inputTable;
    pc_struct.table = table;
elseif(strcmp(pc_struct.type,strcat('BROOKS-COREY')))
    pc_struct.pd = str2double(line(3));
    pc_struct.lambda = str2double(line(4));
elseif(strcmp(pc_struct.type,strcat('BROOKS-COREY-SEPARATE-SWC')))
    pc_struct.Swc = str2double(line(3));
    pc_struct.Sor = str2double(line(4));  
    pc_struct.pd = str2double(line(5));
    pc_struct.lambda = str2double(line(6)); 
elseif(strcmp(pc_struct.type,strcat('VAN-GENUCHTEN')))
    pc_struct.alpha = str2double(line(3));
    pc_struct.n = str2double(line(4));
elseif(strcmp(pc_struct.type,strcat('VAN-GENUCHTEN-SEPARATE-SWC')))
    pc_struct.Swc = str2double(line(3));
    pc_struct.Sor = str2double(line(4));  
    pc_struct.alpha = str2double(line(5));
    pc_struct.n = str2double(line(6)); 
elseif(strcmp(pc_struct.type,strcat('SKJAEVELAND')))
    pc_struct.cwi = str2double(line(3));
    pc_struct.coi = str2double(line(4));
    pc_struct.awi = str2double(line(5));
    pc_struct.aoi = str2double(line(6));
elseif(strcmp(pc_struct.type,strcat('MODIFIED-SKJAEVELAND')))
    pc_struct.cwi = str2double(line(3));
    pc_struct.coi = str2double(line(4));
    pc_struct.ri = str2double(line(5));
    pc_struct.bi = str2double(line(6));
    pc_struct.Swd = str2double(line(7));
    pc_struct.Sod = str2double(line(8));
elseif(strcmp(pc_struct.type,strcat('MODIFIED-SKJAEVELAND-MASALMEH')))
    pc_struct.cwi = str2double(line(3));
    pc_struct.coi = str2double(line(4));
    pc_struct.awi = str2double(line(5));
    pc_struct.aoi = str2double(line(6));
    pc_struct.sw_cutoff = str2double(line(7));
    pc_struct.b = str2double(line(8));
elseif(strcmp(pc_struct.type,strcat('MODIFIED-SKJAEVELAND-MASALMEH-SEPARATE-SWC')))
    pc_struct.Swc = str2double(line(3));
    pc_struct.Sor = str2double(line(4));  
    pc_struct.cwi = str2double(line(5));
    pc_struct.coi = str2double(line(6));
    pc_struct.awi = str2double(line(7));
    pc_struct.aoi = str2double(line(8));
    pc_struct.sw_cutoff = str2double(line(9));
    pc_struct.b = str2double(line(10));
elseif(strcmp(pc_struct.type,strcat('ZERO')))
    pc_struct.zero = true;
elseif(strcmp(pc_struct.type,strcat('LET-IMBIBITION')))
    pc_struct.sw_pc0 = str2double(line(3));
    pc_struct.max_pc = str2double(line(4));
    pc_struct.min_pc = str2double(line(5));
    pc_struct.spontaneous_multiplier = str2double(line(6));
    pc_struct.forced_multiplier = str2double(line(7));
    pc_struct.L_spont = str2double(line(8));
    pc_struct.E_spont = str2double(line(9));
    pc_struct.T_spont = str2double(line(10));
    pc_struct.L_forced = str2double(line(11));
    pc_struct.E_forced = str2double(line(12));
    pc_struct.T_forced = str2double(line(13));                
elseif(strcmp(pc_struct.type,strcat('LET-DRAINAGE')))
    pc_struct.max_pc = str2double(line(3));
    pc_struct.entry_pc = str2double(line(4));
    pc_struct.entry_multiplier = str2double(line(5));
    pc_struct.forced_multiplier = str2double(line(6));
    pc_struct.L_entry = str2double(line(7));
    pc_struct.E_entry = str2double(line(8));
    pc_struct.T_entry = str2double(line(9));
    pc_struct.L_forced = str2double(line(10));
    pc_struct.E_forced = str2double(line(11));
    pc_struct.T_forced = str2double(line(12));
elseif(strcmp(pc_struct.type,strcat('LET-DRAINAGE-SEPARATE-SWC')))
    pc_struct.swc = str2double(line(3));
    pc_struct.max_pc = str2double(line(4));
    pc_struct.entry_pc = str2double(line(5));
    pc_struct.entry_multiplier = str2double(line(6));
    pc_struct.forced_multiplier = str2double(line(7));
    pc_struct.L_entry = str2double(line(8));
    pc_struct.E_entry = str2double(line(9));
    pc_struct.T_entry = str2double(line(10));
    pc_struct.L_forced = str2double(line(11));
    pc_struct.E_forced = str2double(line(12));
    pc_struct.T_forced = str2double(line(13));
end 