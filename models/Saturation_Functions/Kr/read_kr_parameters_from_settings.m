function kr_struct = read_kr_parameters_from_settings(kr_struct, line)

if(strcmp(kr_struct.type,strcat('TABLE')))
    kr_struct.fullFile = line{3};
    [filePath,fileName,fileExt] = fileparts(kr_struct.fullFile);
    kr_struct.filePath = filePath;
    kr_struct.fileName = strcat(fileName,fileExt);
    [~,table] = ImportTable(kr_struct.fullFile);                    
%     [~,inputTable,~,table] = ImportTable2(kr.fullFile);
%     kr.inputTable = inputTable;
    kr_struct.table = table;
elseif(strcmp(kr_struct.type,strcat('MODIFIED-COREY')))                   
    kr_struct.Swc     = str2double(line(3));
    kr_struct.Sor     = str2double(line(4));
    kr_struct.krwSor  = str2double(line(5));
    kr_struct.kroSwc  = str2double(line(6));
    kr_struct.nW      = str2double(line(7));
    kr_struct.nNW     = str2double(line(8)); 
elseif(strcmp(kr_struct.type,strcat('MODIFIED-COREY-MASALMEH')))  
    kr_struct.Swc     = str2double(line(3));
    kr_struct.Sor     = str2double(line(4));
    kr_struct.krwSor  = str2double(line(5));
    kr_struct.kroSwc  = str2double(line(6));
    kr_struct.nW      = str2double(line(7));
    kr_struct.nNW     = str2double(line(8)); 
    kr_struct.cW       = str2double(line(9)); 
    kr_struct.cNW       = str2double(line(9));
elseif(strcmp(kr_struct.type,strcat('BROOKS-COREY')))
    kr_struct.Swc     = str2double(line(3));                    
    kr_struct.Sor     = str2double(line(4));
    kr_struct.krwSor  = str2double(line(5));
    kr_struct.kroSwc  = str2double(line(6));
    kr_struct.lambda  = str2double(line(7));
elseif(strcmp(kr_struct.type,strcat('BURDINE')))
    kr_struct.krwSor  = str2double(line(3));
    kr_struct.kroSwc  = str2double(line(4));   
elseif(strcmp(kr_struct.type,strcat('VAN-GENUCHTEN-BURDINE')))
    kr_struct.Swc     = str2double(line(3));                    
    kr_struct.Sor     = str2double(line(4));
    kr_struct.krwSor  = str2double(line(5));
    kr_struct.kroSwc  = str2double(line(6));
    kr_struct.n  = str2double(line(7));
elseif(strcmp(kr_struct.type,strcat('LET')))
    kr_struct.Swc    = str2double(line(3));
    kr_struct.Sor    = str2double(line(4));                    
    kr_struct.krwSor = str2double(line(5));
    kr_struct.kroSwc = str2double(line(6));
    kr_struct.Lw     = str2double(line(7));
    kr_struct.Lnw    = str2double(line(8));
    kr_struct.Ew     = str2double(line(9));
    kr_struct.Enw    = str2double(line(10));
    kr_struct.Tw     = str2double(line(11));
    kr_struct.Tnw    = str2double(line(12)); 
end