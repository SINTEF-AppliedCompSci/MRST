function a = paramXmlToStruct(file)
xroot=parseXML(file);
a=[];
a=readLevel(xroot, a);
fnames=fields(a);
assert(numel(fnames)==1);
assert(strcmp(fnames{1},'root'));
a=a.root;
end
function a=readLevel(child, a)
for j=1:numel(child)
   atrib={child.Attributes};
   for k=1:numel(atrib)
      if(~isempty(atrib{k}))
         assert(strcmp(atrib{k}(1).Name,'name'))
         if(numel(atrib{k})==3)
            assert(strcmp(atrib{k}(2).Name,'type'))
            assert(strcmp(atrib{k}(3).Name,'value'))
            a.(atrib{k}(1).Value) = atrib{k}(3).Value;
         elseif(numel(atrib{k})==1)            
            assert(strcmp(child.Name,'ParameterGroup'));
            mychild=child.Children;
            if(~isfield(a,atrib{k}(1).Value))
                a.(atrib{k}(1).Value)=[];
            end
            for i=1:numel(mychild)
                a.(atrib{k}(1).Value)=readLevel(mychild(i), a.(atrib{k}(1).Value));
            end
         else
            error('error')
         end
            
      end
   end
end

%mychild=child.Children;
%for i=1:numel(mychild)
%   a=readLevel(mychild(i), a);
%end
end

