function [model_flooding, model_cent] = historymatch(model_flooding)
%
% DESCRIPTION: main history match function for the single objective history
%              matching using fmincon function from MATLAB
%
% SYNOPSIS:
%   [model_flooding, model_cent] = historymatch(model_flooding)
%
% PARAMETERS:
%   model_flooding - struct containing following fields:
%   - history_match - including history matching boundaries and initial
%   points (in case of a simultaneous history matching, relavant info
%   should be included as well)
%   - experiment - heterogeneity can be included here
%
% RETURNS:
%   model_flooding - struct containing following fields:
%   - history_match - includes the history matching results
%   model_cent - secondary struct in case of simultaneous history matching
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
    lb = model_flooding.history_match.lb;
    ub = model_flooding.history_match.ub;
    x0 = model_flooding.history_match.x0;
    multi_point = model_flooding.history_match.multi_point;
    if multi_point, rnd_num = model_flooding.history_match.multipoint_rnd_num; end
    UseParallel = model_flooding.history_match.UseParallel;
    OptimalityTolerance = model_flooding.history_match.OptimalityTolerance;
    StepTolerance = model_flooding.history_match.StepTolerance;
    MaxFunctionEvaluations = model_flooding.history_match.MaxFunctionEvaluations;
    ScaleProblem = model_flooding.history_match.ScaleProblem;
    algorithm = model_flooding.history_match.algorithm;
    % to test the mcmc objective
%     model_flooding.history_match.mcmc = true;
    
    %setup centrifuge model
    obj_fun_type = model_flooding.history_match.obj_fun;
    if strcmp(obj_fun_type,'Simultaneous')
        appDir = model_flooding.history_match.Cent_file_path; 
        name = model_flooding.history_match.Cent_file_name;
        temp_historymatch_variable_hold = model_flooding.history_match;
        cent_directory = fullfile(model_flooding.history_match.Cent_file_path, ...
            model_flooding.history_match.Cent_file_name);
        model_cent = Configure(cent_directory);
        model_cent = CreateGrid(model_cent);
        model_cent = CreateRock(model_cent);
        model_cent = CreatePc(model_cent);
        model_cent = CreateKr(model_cent);
        model_cent.history_match = temp_historymatch_variable_hold;
        % to test the mcmc objective
%         model_cent.history_match.mcmc = true;
    else 
        model_cent = [];
    end
    
    
    % Pass fixed parameters to objfun
    if model_flooding.experiment.rock.heterogeneous
        obj_fun = @(x) objectivefun_sync_heterogenous(x, model_flooding, model_cent);
    else
        obj_fun = @(x) objectivefun_sync(x, model_flooding, model_cent);
    end
    
    % Set up shared variables with OUTFUN
    history.x = [];
    history.fval = [];

    filePath = strcat(pwd);
    fileName = 'historymatch_output.txt';
    fileID = fopen(fullfile(filePath,fileName),'w');
    fprintf(fileID,'%s\n','Algorithm iterations, Objective function iterations, Objective function value,Best input variable values');
    modified_outputfun = @(x,optimValues,state) outputFcn(x,optimValues,state,fileID);          


    % Set nondefault solver options
    options = optimoptions('fmincon','FunValCheck','on','OptimalityTolerance', OptimalityTolerance,...
        'Display','iter','UseParallel',UseParallel, 'FiniteDifferenceType', 'central'...
        ,'PlotFcn', {'optimplotx', 'optimplotfunccount', 'optimplotfval'},...
        'ScaleProblem', ScaleProblem, 'Diagnostics', 'off',...
        'Algorithm',algorithm,'MaxFunctionEvaluations', MaxFunctionEvaluations,...
        'StepTolerance', StepTolerance,'OutputFcn', modified_outputfun);

    % setting up inequalities
    [A, b] = setup_A_b(model_flooding);

    % Solve
    problem = createOptimProblem('fmincon','objective',...
        obj_fun, 'x0', x0, 'Aineq', A, 'bineq', b, 'lb', lb, 'ub', ub,'options',options);

    ms = MultiStart('Display','iter', 'StartPointsToRun','bounds-ineqs',...
        'PlotFcn', {'gsplotbestf','gsplotfunccount'});
    
    if multi_point
        [x,fval] = run(ms,problem,rnd_num + 1);
    else
        [x,fval,~,~,~,~,hessian] = fmincon(obj_fun,x0,A,b,[],[],lb,ub,[],options);
        model_flooding.history_match.hessian = hessian;
        model_flooding.history_match.error = sqrt(abs(diag(inv(hessian))));
    end

    model_flooding.history_match.x = x;
    model_flooding.history_match.fval = fval;
    model_flooding.history_match.history = history;

    % Clear variables
    clearvars fun options
    fclose(fileID);
    
    
    % writing the hm results in the input excel file
    hm_template_name = model_flooding.history_match.hm_template_name;
    hm_template_path = model_flooding.history_match.hm_template_path;
    writematrix('History_match_results',fullfile(hm_template_path,hm_template_name)...
        ,'Sheet',1,'Range','M1')    
    writematrix(model_flooding.history_match.x',fullfile(hm_template_path,hm_template_name)...
        ,'Sheet',1,'Range','M2')
    
    function stop = outputFcn(x,optimValues,state,fileID)
         stop = false;
         switch state
             case 'init'
             case 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
                 history.fval = [history.fval; optimValues.fval];
                 history.x = [history.x; x];
                 if(~isempty(fileID))
                    fprintf(fileID,'%12.6f,',optimValues.iteration,optimValues.funccount, optimValues.fval, x);
                    fprintf(fileID,'\n');
                 end
             case 'done'
             otherwise
         end
    end
end