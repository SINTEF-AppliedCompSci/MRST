function state = initGasPhasePressures(model, state, partialPressures, nI);


givenNames = partialPressures(1:2:end);
givenValues = partialPressures(2:2:end);

assert(numel(givenNames) == model.nG, 'Number of gas phase pressures given does not match number of gas phases in the system.');

ind = zeros(1, model.nG);
for i = 1 : model.nG
    ind = ind + strcmpi(model.gasNames{i}, givenNames);
end
ind = logical(ind);

missingNames = model.gasNames(~ind);
nMN = numel(missingNames);

if nMN > 0
    iwant = '';
    for i = 1 : nMN
        if i == nMN
            iwant = [iwant missingNames{:}];
        else
           iwant = [iwant , ' , '];
        end
    end
    
    error( ['The partial pressure of ' iwant ' are missing from partialPressures.'] );
end

state.partialPressures = zeros(1, model.nS);

for i = 1 : model.nG
    ind = strcmpi(model.gasNames{i}, givenNames);
    state.partialPressures(i) = givenValues{ind};
end

if size(state.partialPressures,1) == 1
    state.partialPressures = repmat(state.partialPressures,nI,1);
elseif size(state.partialPressures,1) ~= nI
    error('Size of partialPressures is not one and does not match the size of userInput.');
end


end

