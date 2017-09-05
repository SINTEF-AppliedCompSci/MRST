function chemicalSystemPrintFunction(model, varargin )

            fprintf('\n');
            
            l1 = max(cellfun(@length, model.rxns))+1;
            l2 = max(cellfun(@length, model.AllCompNames))+1;

            sep = repmat('=', 1, l2*numel(model.AllCompNames) + l1 + numel(model.AllCompNames) + 2);
            fprintf('%s\n', sep);
            
            % get the names right
            for i = 1 : model.nMC
            	mcnames{i} = ['sum(' model.MasterCompNames{i} ')'];
            end
            for i = 1 : model.nR
                mcnames{model.nMC+i} = model.rxns{i};
            end

            % make column header
            compnamestr =['|%-' num2str(l1) 's|'];
            for i = 1 : numel(model.AllCompNames)
                compnamestr = [compnamestr '%' num2str(l2) 's|'];
            end
            fprintf([compnamestr '\n'], 'Equations', model.AllCompNames{:});

            fprintf('%s\n', sep);
            % make the rows of the print
            
            massnamestr ='';
            for i = 1 : model.nMC + model.nR
                for j = 1 : numel(model.AllCompNames) + 1
                    if j == 1
                        fl = ['|%-' num2str(l1) 's'];
                    else
                        fl = ['|%' num2str(l2) 'g'];
                    end

                    massnamestr = [massnamestr '' fl ''];
                end
                if i <= model.nMC
                    in = mat2cell(model.AllContributionMatrix(i,:), 1, ones(1,model.nC+model.nG+model.nS) );
                else
                    in = mat2cell(model.AllReactionMatrix(i-model.nMC,:), 1, ones(1,model.nC+model.nG+model.nS) );
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
                for j = 1 : numel(model.AllCompNames) + 1
                    if j == 1
                        fl = ['|%-' num2str(l1) 's'];
                    else
                        fl = ['|%' num2str(l2) 'g'];
                    end

                    massnamestr = [massnamestr '' fl ''];
                end
                in = mat2cell(model.AllCombinationMatrix(i,:), 1, ones(1,model.nC+model.nG+model.nS) );
                fprintf([massnamestr, '|\n'], model.CombinationNames{i}, in{:});
                massnamestr ='';
            end
            if model.nLC >0
                fprintf('%s\n', sep);
            end
            
            fprintf('\n');


end

