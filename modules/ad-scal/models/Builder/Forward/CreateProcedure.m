function procedure = CreateProcedure(table, processType)
%
% DESCRIPTION: creates procedure table from the schedule table
%              which is used to run the simulation
%
% SYNOPSIS:
%   procedure = CreateProcedure(table, processType)
%
% PARAMETERS:
%   table - schedule table read from a .txt file
%   processType - can be "SS", "USS" or "CENT"
%
% RETURNS:
%   procedure - table used in the Run.m module for simulation
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
if(strcmpi(processType,"SS") || strcmpi(processType,"USS"))
    headers = table.Properties.VariableNames;
    [headers,~] = SplitHeaders(headers);

    timeFrom  = table{1:end-1,1};
    timeTo    = table{2:end,1}; 
    totalRate = table{1:end-1,2};
    oilFrac   = table{1:end-1,3};
    watFrac   = table{1:end-1,4}; 
    oilRate   = table{1:end-1,2} .* table{1:end-1,3}; 
    watRate   = table{1:end-1,2} .* table{1:end-1,4};
    oilCum    = cumsum(oilRate .* (timeTo - timeFrom));
    watCum    = cumsum(watRate .* (timeTo - timeFrom)); 

    data = [timeFrom, timeTo, totalRate, oilFrac, watFrac, ...
                      oilRate, watRate, oilCum, watCum];

    % table headers
    oilRateHeader  = 'oil rate'; oilRateQuantity = 'rate';
    watRateHeader  = 'water rate'; watRateQuantity = 'rate';
    oilCumHeader   = 'oil cum'; oilCumQuantity = 'volume';
    watCumHeader   = 'water cum'; watCumQuantity = 'volume';   
    newHeaders{1}  = strcat(headers{1},' [',Unit(headers{1}),']',': from');
    newHeaders{2}  = strcat(headers{1},' [',Unit(headers{1}),']',': to');
    newHeaders{3}  = strcat(headers{2},' [',Unit(headers{2}),']');
    newHeaders{4}  = strcat(headers{3},' [',Unit(headers{3}),']');
    newHeaders{5}  = strcat(headers{4},' [',Unit(headers{4}),']');
    newHeaders{6}  = strcat(oilRateHeader,' [',Unit(oilRateQuantity),']');
    newHeaders{7}  = strcat(watRateHeader,' [',Unit(watRateQuantity),']');
    newHeaders{8}  = strcat(oilCumHeader,' [',Unit(oilCumQuantity),']');
    newHeaders{9}  = strcat(watCumHeader,' [',Unit(watCumQuantity),']');  
end
if(strcmpi(processType,"CENT") || strcmpi(processType,"Centrifuge"))
    headers = table.Properties.VariableNames;
    [headers,~] = SplitHeaders(headers);

    timeFrom  = table{1:end-1,1};
    timeTo    = table{2:end,1}; 
    acceleration = table{1:end-1,2}; 

    data = [timeFrom, timeTo, acceleration];

    % table headers 
    newHeaders{1}  = strcat(headers{1},' [',Unit(headers{1}),']',': from');
    newHeaders{2}  = strcat(headers{1},' [',Unit(headers{1}),']',': to');
    newHeaders{3}  = strcat(headers{2},' [',Unit(headers{2}),']');
end
procedure = array2table(data,'VariableNames',newHeaders);