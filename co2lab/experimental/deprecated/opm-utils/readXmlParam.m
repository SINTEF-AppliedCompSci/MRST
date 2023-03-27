function a = readXmlParam(file)
xroot=parseXML(file);
a=[];
a=readLevel(xroot, a);
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
            a.(atrib{k}(1).Value)='ParamGroup';
         else
            error('error')
         end
            
      end
   end
end
mychild=child.Children;
for i=1:numel(mychild)
   a=readLevel(mychild(i), a);
end
end

