function [eqs, names, types] = equationsCompositionGuess(comps, masterComps, model)
    
%             comps = cellfun(@(x) x*litre/mol, comps,'UniformOutput', false);
%             masterComps = cellfun(@(x) x*litre/mol, masterComps,'UniformOutput', false);
%             
% begin solving equation

    CNames = model.CompNames;

    compInd = cellfun(@(x) isempty(x), regexpi(CNames, 'psi'));
    CNames = CNames(compInd);
    nC = numel(CNames);

    CM = model.CompositionMatrix;
    CM(:,~compInd) = [];
            
    eqs   = cell(1, model.nMC);
    names = cell(1, model.nMC);
    types = cell(1, model.nMC);


    for i = 1 : model.nMC
        eqs{i} = - masterComps{i};
        for k = 1 : nC
            eqs{i} = eqs{i} + CM(i,k).*comps{k};
        end
        names{i} = ['Conservation of ', model.MasterCompNames{i}] ;
    end

    [types{:}] = deal('cell');
    
end

