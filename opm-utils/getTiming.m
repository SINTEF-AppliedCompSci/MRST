function myvec= getTiming(param,myfield,mytype)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
param{1}=convertTiming(param{1});
myvec=zeros(numel(param),numel(param{1}.timing.(myfield).(mytype)));
for i=1:numel(param)
   param{i}=convertTiming(param{i});
   myvec(i,:)=param{i}.timing.(myfield).(mytype)(:);
end
end

