function chemicalSystemPrintFunction(model, varargin )

            fprintf('\n');
            
            l1 = max(cellfun(@length, model.rxns))+1;
            l2 = max(cellfun(@length, model.CompNames))+1;

            sep = repmat('=', 1, l2*model.nC + l1 + model.nC + 2);
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
            for i = 1 : model.nC
                compnamestr = [compnamestr '%' num2str(l2) 's|'];
            end
            fprintf([compnamestr '\n'], 'Equations', model.CompNames{:});

            fprintf('%s\n', sep);
            % make the rows of the print
            
            massnamestr ='';
            for i = 1 : model.nMC + model.nR
                for j = 1 : model.nC + 1
                    if j == 1
                        fl = ['|%-' num2str(l1) 's'];
                    else
                        fl = ['|%' num2str(l2) 'g'];
                    end

                    massnamestr = [massnamestr '' fl ''];
                end
                if i <= model.nMC
                    in = mat2cell(model.CompositionMatrix(i,:), 1, ones(1,model.nC) );
                else
                    in = mat2cell(model.ReactionMatrix(i-model.nMC,:), 1, ones(1,model.nC) );
                end
                if i == model.nMC+1
                    fprintf('%s\n', sep);
                end
                fprintf([massnamestr, '|\n'], mcnames{i}, in{:});
                massnamestr ='';
            end
            fprintf('%s\n', sep);
            fprintf('\n');


end

