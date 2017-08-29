function model = initSolidPhaseDensities(model, solidDensities);


givenNames = solidDensities(1:2:end);
givenValues = solidDensities(2:2:end);

assert(numel(givenNames) == model.nS, 'Number of solid phase densities given does not match number of solid phases in the system.');

ind = zeros(1, model.nS);
for i = 1 : model.nS
    ind = ind + strcmpi(model.SolidNames{i}, givenNames);
end
ind = logical(ind);

missingNames = model.SolidNames(~ind);
nMN = numel(missingNames);

if nMN > 0
    iwant = '';
    for i = 1 : nMN
        if i == nMN
            iwant = [iwant missingNames];
        else
           iwant = [iwant , ' , '];
        end
    end
    
    error( ['The solid densities of ' iwant ' are missing from solidDensities.'] );
end

model.solidDensities = zeros(1, model.nS);

for i = 1 : model.nS
    ind = strcmpi(model.SolidNames{i}, givenNames);
    model.solidDensities(i) = givenValues{ind};
end

end

