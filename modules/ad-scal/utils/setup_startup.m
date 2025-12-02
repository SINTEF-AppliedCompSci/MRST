function [procedure, schedule] = setup_startup(schedule, procedure, value_RPM, unit_RPM)

schedule.startupRPM.inputValue = value_RPM;
schedule.startupRPM.inputUnit = string(unit_RPM);  
schedule.startupRPM.value = (value_RPM * Convert(unit_RPM))...
    ^ 2 * schedule.centRad.value;
headers = procedure.Properties.VariableNames;
procedure_array = table2array(procedure); 
% schedule_in_rpm
SIR = sqrt(procedure_array(:,3)/schedule.centRad.value)...
    / Convert('rpm');
% procedure_table_in_rpm
PTIR = [procedure_array(:,1:2),SIR];
% start_up_rpm
SUR = sqrt(schedule.startupRPM.value/schedule.centRad.value)...
    / Convert('rpm');
% startupSlope
stSL = SUR / schedule.startupPeriod.value;
% added_schedule_row                 
ASRow = [];
% interval_between_steps
IBS = 20; % seconds
for i = 1 : height(PTIR)
    if i == 1
        % time_needed
        TN = PTIR(i,3) / stSL;
    else
        TN = ( PTIR(i,3) - PTIR(i-1,3)) / stSL; 
    end
    % residual_value
    RV = rem(TN, IBS);
    % number_of_steps
    NOS = floor(TN / IBS);
    % added_schedule_times
    AST = linspace(PTIR(i,1),PTIR(i,1) + NOS * IBS, NOS + 1);
    if i == 1
        % added_schedule_rpm
        ASrpm = linspace(0 ,NOS * IBS, NOS + 1) * stSL;
    else
        ASrpm = linspace(0 ,NOS * IBS, NOS + 1) * stSL + ...
            PTIR(i-1,3);   
    end
    for j = 1 : length(AST) - 1
        ASRow = [ASRow; ...
            [AST(j), AST(j + 1), ASrpm(j + 1)]];
    end
    if not(RV == 0)
        ASRow = [ASRow; ...
            [AST(end), AST(end) + RV, ASrpm(end) + RV * stSL]];
    end
    ASRow = [ASRow; ...
        [AST(end) + RV, PTIR(i,2), PTIR(i,3)]];
end   
procedure = [ASRow(:,1:2),(ASRow(:,3)*Convert('rpm')) .^2 .* ...
    schedule.centRad.value];               
procedure = array2table(procedure,'VariableNames',headers);