function kr_struct = read_kr_parameters_from_settings(kr_struct, splitted_lin, path_to_settings)

if(strcmp(kr_struct.type,strcat('TABLE')))
    kr_struct.fullFile = abs_path(splitted_lin{3}, path_to_settings);
    [filePath,fileName,fileExt] = fileparts(kr_struct.fullFile);
    kr_struct.filePath = filePath;
    kr_struct.fileName = strcat(fileName,fileExt);
    [~,table] = ImportTable(kr_struct.fullFile);                    
%     [~,inputTable,~,table] = ImportTable2(kr.fullFile);
%     kr.inputTable = inputTable;
    kr_struct.table = table;
elseif(strcmp(kr_struct.type,strcat('MODIFIED-COREY')))                   
    kr_struct.Swc     = str2double(splitted_lin{3});
    kr_struct.Sor     = str2double(splitted_lin{4});
    kr_struct.krwSor  = str2double(splitted_lin{5});
    kr_struct.kroSwc  = str2double(splitted_lin{6});
    kr_struct.nW      = str2double(splitted_lin{7});
    kr_struct.nNW     = str2double(splitted_lin{8}); 
elseif(strcmp(kr_struct.type,strcat('MODIFIED-COREY-MASALMEH')))  
    kr_struct.Swc     = str2double(splitted_lin{3});
    kr_struct.Sor     = str2double(splitted_lin{4});
    kr_struct.krwSor  = str2double(splitted_lin{5});
    kr_struct.kroSwc  = str2double(splitted_lin{6});
    kr_struct.nW      = str2double(splitted_lin{7});
    kr_struct.nNW     = str2double(splitted_lin{8}); 
    kr_struct.cW      = str2double(splitted_lin{9}); 
    kr_struct.cNW     = str2double(splitted_lin{10});
elseif(strcmp(kr_struct.type,strcat('BROOKS-COREY')))
    kr_struct.Swc     = str2double(splitted_lin{3});                    
    kr_struct.Sor     = str2double(splitted_lin{4});
    kr_struct.krwSor  = str2double(splitted_lin{5});
    kr_struct.kroSwc  = str2double(splitted_lin{6});
    kr_struct.lambda  = str2double(splitted_lin{7});
elseif(strcmp(kr_struct.type,strcat('BURDINE')))
    kr_struct.krwSor  = str2double(splitted_lin{3});
    kr_struct.kroSwc  = str2double(splitted_lin{4});   
elseif(strcmp(kr_struct.type,strcat('VAN-GENUCHTEN-BURDINE')))
    kr_struct.Swc     = str2double(splitted_lin{3});                    
    kr_struct.Sor     = str2double(splitted_lin{4});
    kr_struct.krwSor  = str2double(splitted_lin{5});
    kr_struct.kroSwc  = str2double(splitted_lin{6});
    kr_struct.n       = str2double(splitted_lin{7});
elseif(strcmp(kr_struct.type,strcat('LET')))
    kr_struct.Swc    = str2double(splitted_lin{3});
    kr_struct.Sor    = str2double(splitted_lin{4});                    
    kr_struct.krwSor = str2double(splitted_lin{5});
    kr_struct.kroSwc = str2double(splitted_lin{6});
    kr_struct.Lw     = str2double(splitted_lin{7});
    kr_struct.Lnw    = str2double(splitted_lin{8});
    kr_struct.Ew     = str2double(splitted_lin{9});
    kr_struct.Enw    = str2double(splitted_lin{10});
    kr_struct.Tw     = str2double(splitted_lin{11});
    kr_struct.Tnw    = str2double(splitted_lin{12}); 
end

end

%% Translate relative pathname to absolute pathname.
function path = abs_path(inc_fn, setting_path)
inc_fn(or(inc_fn == '/', inc_fn == '\')) = filesep;
    if inc_fn(1) ~= filesep
       % Translate relative pathname to absolute pathname.
       path = fullfile(fileparts(setting_path), inc_fn);
    end
end