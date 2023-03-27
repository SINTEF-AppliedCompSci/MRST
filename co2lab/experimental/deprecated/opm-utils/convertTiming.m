function param = convertTiming(param)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
myfields=fieldnames(param);
for i=1:numel(myfields)
   myfield=myfields{i};
   if(iscell(param.(myfield)))
      a=cellfun(@(x) str2num(x),param.(myfield),'UniformOut',false);
      param.(myfield)=[a{:}];
   elseif(ischar(param.(myfield)))
      param.(myfield)=str2num(param.(myfield));
   elseif(isstruct(param.(myfield)))  
      param.(myfield)=convertTiming(param.(myfield));
   end
end
end

