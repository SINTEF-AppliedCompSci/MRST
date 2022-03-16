function sampleStructure = postProcessRockSamples(mainDirectory, ranges, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('saveFigFolder', '', ...
                 'trueRock', struct());
    
    opt = merge_options(opt, varargin{:});
             
    doSave = true;
    if strcmp(opt.saveFigFolder, '')
        doSave = false;
    end
    sampleStructure = getSamplesStructure(mainDirectory);
    
    % Compute mean and variance for each parameter and sample
    
    for it = 1:numel(sampleStructure)
        for subit = 1:numel(sampleStructure{it})
            sampleStructure{it}{subit} = getUnscaledSampleVectors(sampleStructure{it}{subit});
        end
    end
    
    titles = {'Permeabilities', 'Porosities'};
    xlabels = {'Grid cell', 'Grid cell'};
    parameter_names = {'permeabilities', 'porosities'};
    fileExt = {'.fig', '.png', '.eps'}; %, '.pdf'};
    formats = {'fig', 'png', 'epsc'}; %, 'pdf'};
    
    trueFields = {'perm', 'poro'};
    
    for i = 1:numel(ranges)
        fig = figure();
        fig.Position = [100, 350, 1000, 250];
        b = boxchart(sampleStructure{1}{1}.unscaledSampleVectors(ranges{i}, :)');
        hold on
        if numel(sampleStructure) > 1
            b = [b; boxchart(sampleStructure{2}{1}.unscaledSampleVectors(ranges{i}, :)')];
        end
        xlabel(xlabels{i});
        title(titles{i});
        if i == 1
            fig.CurrentAxes.YAxis.Scale = 'log';
        end
                
        if isfield(opt.trueRock, trueFields{i})
            trueVals = getfield(opt.trueRock, trueFields{i});
            b = [b; plot(trueVals, 'xk')];
            if numel(sampleStructure) > 1
                legend(b, 'prior', 'posterior', 'truth');
            else
                legend(b, 'prior', 'truth'); 
            end
        else
            if numel(sampleStructure) > 1
                legend(b, 'prior', 'posterior');
            else
                legend(b, 'prior');
            end
        end
        
        hold off
        if doSave
            for f = 1:numel(fileExt)
                savefilename = fullfile(opt.saveFigFolder, strcat(parameter_names{i}, fileExt{f}))
                saveas(fig, savefilename, formats{f});
            end
        end
    end
    
    saveSample(mainDirectory, 'priorSamples', ...
        sampleStructure{1}{1}.unscaledSampleVectors, ranges);
    
    for i = 2:numel(sampleStructure{1})
        saveSample(mainDirectory, strcat('interimSamples_PostIteration', num2str(i-1)), ...
            sampleStructure{1}{i}.unscaledSampleVectors, ranges);
    end
    
    if numel(sampleStructure) > 1
        saveSample(mainDirectory, 'posteriorSamples', ...
            sampleStructure{2}{1}.unscaledSampleVectors, ranges);
    end    
    
    
    
    
end

function saveSample(mainDirectory, filename, unscaledSamples, ranges)
    fullFilename = fullfile(mainDirectory, filename);
    permeabilites = unscaledSamples(ranges{1}, :);
    porosities = unscaledSamples(ranges{2}, :);
    
    save(fullFilename, 'permeabilites', 'porosities');
end


function sampleStruct = getUnscaledSampleVectors(sampleStruct)
    
    % TODO: Check for CompositeSample and recursively iterate down the
    % rabbit hole
    %for i = 1:numel(sampleStruct.samples.parentSamples)
    %    originalTransformSampleVectors(i) = sampleStruct.samples.parentSamples{i}.transformSampleVectors;
    %    sampleStruct.samples.parentSamples{i}.transformSampleVectors = false;
    %end
    %
    %sampleStruct.unscaledSampleVectors = sampleStruct.samples.getSampleVectors();
    % 
    %for i = 1:numel(sampleStruct.samples.parentSamples)
    %    sampleStruct.samples.parentSamples{i}.transformSampleVectors = originalTransformSampleVectors(i);
    %end
    
    originalTransformSampleVectors = sampleStruct.samples.transformSampleVectors;
    sampleStruct.samples.transformSampleVectors = false;
    sampleStruct.unscaledSampleVectors = sampleStruct.samples.getSampleVectors();
    sampleStruct.samples.transformSampleVectors = originalTransformSampleVectors;
    
end





function sampleStructure = getSamplesStructure(mainDirectory)

    sampleStructure = {};
    
    content = dir(mainDirectory);
    iterationDirs = content([content.isdir]);
    
    itNum = 1;
    for i = 1:numel(iterationDirs)
        if strcmp(iterationDirs(i).name, num2str(itNum))
            
            subitNum = 1;
            subFolder = fullfile(mainDirectory, iterationDirs(i).name);
            subContent = dir(subFolder);
            subDirs = subContent([subContent.isdir]);
            
            subSampleStructure = {};
            
            for j = 1:numel(subDirs)
                if strcmp(subDirs(j).name, num2str(subitNum))
                    sampleCandidatePath = fullfile(subFolder, subDirs(j).name, 'samples.mat');
                    if exist(sampleCandidatePath, 'file')
                        subSampleStructure{subitNum} = load(sampleCandidatePath);
                        subitNum = subitNum + 1;
                    end
                end 
            end
            
            if numel(subSampleStructure) >  0
                sampleStructure{itNum} = subSampleStructure;
                itNum = itNum + 1;
            end
            
        end
    end


end
