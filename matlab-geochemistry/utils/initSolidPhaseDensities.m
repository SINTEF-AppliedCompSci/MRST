function state = initSolidPhaseDensities(model, state, solidDensities,nI)


givenNames = solidDensities(1:2:end);
givenValues = solidDensities(2:2:end);

assert(numel(givenNames) == model.nS, 'Number of solid phase densities given does not match number of solid phases in the system.');

ind = zeros(1, model.nS);
for i = 1 : model.nS
    ind = ind + strcmpi(model.solidNames{i}, givenNames);
end
ind = logical(ind);

missingNames = model.solidNames(~ind);
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
    
    error( ['The solid densities of ' iwant{:} ' are missing from solidDensities.'] );
end

state.solidDensities = zeros(1, model.nS);

for i = 1 : model.nS
    ind = strcmpi(model.solidNames{i}, givenNames);
    state.solidDensities(i) = givenValues{ind};
end

if size(state.solidDensities,1) == 1
    state.solidDensities = repmat(state.solidDensities,nI,1);
elseif size(state.solidDensities,1) ~= nI
    error('Size of solidDensities is not one and does not match the size of userInput.');
end

end

