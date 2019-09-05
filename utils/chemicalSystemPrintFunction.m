function chemicalSystemPrintFunction(model, varargin )

            fprintf('\n');
            
            l1 = max(cellfun(@length, model.rxns))+1;
            l2 = max(cellfun(@length, model.allComponentNames))+1;

            sep = repmat('=', 1, l2*numel(model.allComponentNames) + l1 + numel(model.allComponentNames) + 2);
            fprintf('%s\n', sep);
            
            % get the names right
            for i = 1 : model.nMC
            	mcnames{i} = ['sum(' model.elementNames{i} ')'];
            end
            for i = 1 : model.nR
                mcnames{model.nMC+i} = model.rxns{i};
            end

            % make column header
            compnamestr =['|%-' num2str(l1) 's|'];
            for i = 1 : numel(model.allComponentNames)
                compnamestr = [compnamestr '%' num2str(l2) 's|'];
            end
            fprintf([compnamestr '\n'], 'Equations', model.allComponentNames{:});

            fprintf('%s\n', sep);
            % make the rows of the print
            
            massnamestr ='';
            for i = 1 : model.nMC + model.nR
                for j = 1 : numel(model.allComponentNames) + 1
                    if j == 1
                        fl = ['|%-' num2str(l1) 's'];
                    else
                        fl = ['|%' num2str(l2) 'g'];
                    end

                    massnamestr = [massnamestr '' fl ''];
                end
                if i <= model.nMC
                    in = mat2cell(model.allContributionMatrix(i,:), 1, ones(1,model.nC+model.nG+model.nS+model.nP) );
                else
                    in = mat2cell(model.allReactionMatrix(i-model.nMC,:), 1, ones(1,model.nC+model.nG+model.nS+model.nP) );
                end
                if i == model.nMC+1
                    fprintf('%s\n', sep);
                end
                fprintf([massnamestr, '|\n'], mcnames{i}, in{:});
                massnamestr ='';
            end
            fprintf('%s\n', sep);
            
            % linear combinations
            massnamestr ='';
            for i = 1 : model.nLC
                for j = 1 : numel(model.allComponentNames) + 1
                    if j == 1
                        fl = ['|%-' num2str(l1) 's'];
                    else
                        fl = ['|%' num2str(l2) 'g'];
                    end

                    massnamestr = [massnamestr '' fl ''];
                end
                in = mat2cell([model.allCombinationMatrix(i,:), zeros(1,model.nG+model.nS+model.nP)], 1, ones(1,model.nC+model.nG+model.nS+model.nP) );
                fprintf([massnamestr, '|\n'], model.combinationNames{i}, in{:});
                massnamestr ='';
            end
            if model.nLC >0
                fprintf('%s\n', sep);
            end
            
            fprintf('\n');


end

