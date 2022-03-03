function model_flooding = historymatch_multi_obj(model_flooding)
% <keywords>
%
% Purpose : history match the saturation functions using multi objective
% optimization
%
% Syntax : 
%   model_flooding = historymatch_multi_obj(model_flooding)
%
% Input Parameters :
%   model_flooding: struct containing history matching parameters ( in case
%   of a simultaneous optimization, the struct of the second struct is
%   created based on the specs of this struct
%
% Return Parameters :
%   model_flooding: struct containing the history matching resutls
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
    lb = model_flooding.history_match.lb;
    ub = model_flooding.history_match.ub;
    x0 = model_flooding.history_match.x0;
    UseParallel = model_flooding.history_match.UseParallel;
    algorithm = model_flooding.history_match.algorithm;
    
    %setup centrifuge model
    obj_fun_type = model_flooding.history_match.obj_fun;
    if strcmp(obj_fun_type,'Simultaneous')
        appDir = model_flooding.history_match.Cent_file_path; 
        name = model_flooding.history_match.Cent_file_name;
        temp_historymatch_variable_hold = model_flooding.history_match;
        model_cent = Configure(appDir,name);
        model_cent = CreateGrid(model_cent);
        model_cent = CreateRock(model_cent);
        model_cent = CreatePc(model_cent);
        model_cent = CreateKr(model_cent);
        model_cent.history_match = temp_historymatch_variable_hold;
    else 
        model_cent = [];
    end
    
    % Pass fixed parameters to objfun
    if model_flooding.experiment.rock.heterogeneous
        obj_fun = @(x) objectivefun_sync_multi_obj_heterogeneous(x, model_flooding, model_cent);
    else
        obj_fun = @(x) objectivefun_sync_multi_obj(x, model_flooding, model_cent);
    end

    % Set nondefault solver options
%     options = optimoptions( 'paretosearch' , 'UseParallel',UseParallel, ... 
%         'PlotFcn' ,'psplotparetof','Display','iter');
    obs_no = length(fieldnames(model_flooding.experiment.observation));
    if isfield(model_cent, 'experiment')
       obs_no = obs_no + length(fieldnames(model_cent.experiment.observation));
    end
    objectivesToPlot = 1:obs_no;
    plotfn = @(options,state,flag)gaplotpareto(options,state,flag,objectivesToPlot);
    options = optimoptions('gamultiobj','PlotFcn',{plotfn}, ...
        'UseParallel',UseParallel,'Display','iter','InitialPopulationMatrix', x0);

    % setting up inequalities
    [A, b] = setup_A_b(model_flooding);

    % Solve
    [x_list,fval_list] = gamultiobj(obj_fun,numel(x0),A,b,[],[],lb,ub,options);
%     [x,fval] = paretosearch(obj_fun,numel(x0),A,b,[],[],lb,ub,[],options);

    model_flooding.history_match.x_list = x_list;
    model_flooding.history_match.fval_list = fval_list;
end