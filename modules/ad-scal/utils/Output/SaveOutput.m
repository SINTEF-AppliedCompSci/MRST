function SaveOutput(model,fullFile)   
%
% DESCRIPTION: saves the simulation results into an output file
%
% SYNOPSIS:
%   SaveOutput(model, fullFile) 
%
% PARAMETERS:
%   model - struct containing following fields:
%   - dynamic: information of the simulation timespteps
%   - output: settings related to the output reports
%   fullfile - string to the output file destination
%
% RETURNS:
%   outputs the simulation results in a report .txt and excel file
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
params = model.dynamic.params;
output = model.output;

% display units
displayUnits  = DisplayUnits(model);    
displaySat    = displayUnits.displaySat;
displayTime   = displayUnits.displayTime;
displayLength = displayUnits.displayLength;
displayPress  = displayUnits.displayPress;
displayRate   = displayUnits.displayRate;
displayVolume = displayUnits.displayVolume;

% unit converter
satFac   = Convert(displaySat);
timeFac  = Convert(displayTime);
pressFac = Convert(displayPress);
lenFac   = Convert(displayLength);
rateFrac = Convert(displayRate);
volFac   = Convert(displayVolume);
    
time        = params.cumScheduleSteps ./ timeFac;
SwAvg       = params.SwAvg ./ satFac;
pDiff       = params.pDiff ./ pressFac;
qw_inj      = params.qinj(:,1)  ./ rateFrac;
qo_inj      = params.qinj(:,2)  ./ rateFrac;
Qw_inj      = params.Qinj(:,1)  ./ volFac;
Qo_inj      = params.Qinj(:,2)  ./ volFac;
qw_prod     = params.qprod(:,1) ./ rateFrac;
qo_prod     = params.qprod(:,2) ./ rateFrac;
Qw_prod     = params.Qprod(:,1) ./ volFac;
Qo_prod     = params.Qprod(:,2) ./ volFac;
qw_prod_net = params.qp_net(:,1) ./ rateFrac;
qo_prod_net = params.qp_net(:,2) ./ rateFrac;
Qw_prod_net = params.Qp_net(:,1) ./ volFac;
Qo_prod_net = params.Qp_net(:,2) ./ volFac;
PVI         = params.PVI;      

colWidth = 12;
numDigits = 4;
contentFormat = strcat('%',num2str(colWidth),'.',num2str(numDigits),'f');
headerFormat = strcat('%',num2str(colWidth),'s');

headers = []; units = [];
if(isfield(output,'time'))
    if(output.time.include)
    headers = [headers,{'Time'}];
    units = [units,{displayTime}];
    end
end
if(isfield(output,'swavg'))
    if(output.swavg.include)
        headers = [headers,{'SwAvg'}];
        units = [units,{displaySat}];
    end
end
if(isfield(output,'deltaP'))
    if(output.deltaP.include)
    headers = [headers,{'pDiff'}];
    units = [units,{displayPress}];
    end
end
if(isfield(output,'inj'))
    if(output.inj.include)
    headers = [headers,{'qw_inj'},{'qo_inj'},{'Qw_inj'},{'Qo_inj'}]; 
    units = [units,{displayRate},{displayRate},{displayVolume},{displayVolume}];
    end
end
if(isfield(output,'prod'))
    if(output.prod.include)            
    headers = [headers,{'qw_prod'},{'qo_prod'},{'Qw_prod'},{'Qo_prod'},{'Qw_prod_net'},{'Qo_prod_net'}]; 
    units = [units,{displayRate},{displayRate},{displayVolume},{displayVolume},{displayVolume},{displayVolume}];                  
    end
end
headers = string(headers);        
units = string(units);    

file  = fopen(fullFile,'w');
idy   = length(headers) * colWidth + 5;
symbol = '-';
line1 = repelem(symbol,idy);
fprintf(file,'%s\n',line1);

fprintf(file,headerFormat,headers);
fprintf(file,'\n');
fprintf(file,headerFormat,units);
fprintf(file,'\n');

symbol = '-';
line2 = repelem(symbol,idy);
fprintf(file,'%s\n',line2);

for i = 1 : length(time)
    if(isfield(output,'time'))
        if(output.time.include), fprintf(file,contentFormat,time(i)); end
    end
    if(isfield(output,'swavg'))
        if(output.swavg.include), fprintf(file,contentFormat,SwAvg(i)); end
    end
    if(isfield(output,'deltaP'))
        if(output.deltaP.include), fprintf(file,contentFormat,pDiff(i)); end
    end
    if(isfield(output,'inj'))
        if(output.inj.include)
            fprintf(file,contentFormat,qw_inj(i));
            fprintf(file,contentFormat,qo_inj(i));
            fprintf(file,contentFormat,Qw_inj(i));
            fprintf(file,contentFormat,Qo_inj(i));
        end
    end
    if(isfield(output,'prod'))
        if(output.prod.include)
            fprintf(file,contentFormat,qw_prod(i));
            fprintf(file,contentFormat,qo_prod(i));
            fprintf(file,contentFormat,Qw_prod(i));
            fprintf(file,contentFormat,Qo_prod(i));
            fprintf(file,contentFormat,Qw_prod_net(i));
            fprintf(file,contentFormat,Qo_prod_net(i));            
        end
    end
    fprintf(file,'\n');           
end
fclose(file);        