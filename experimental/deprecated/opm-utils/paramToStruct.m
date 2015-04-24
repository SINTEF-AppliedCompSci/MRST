function s=paramToStruct(myfile)
% last line have to have lineshift
myfile=fileread(myfile);
k = regexp(myfile,'([./\w]+)=([-./\w]+)\n','tokens');
k = vertcat(k{:}) .';
s=[];
for i=1:size(k,2)
   b=k{1,i};
   ind=find(b=='/');
   if(isempty(ind))
      ind=[0,numel(b)+1]; 
   else
    if(ind(1)~=1)
          ind=[0,ind];%numel(b)+1]; 
    end
    ind=[ind,numel(b)+1];
   end
    
   names={};
   for jj=1:numel(ind)-1
     names{jj}=b((ind(jj)+1):(ind(jj+1)-1));         
      %s.(name)=[];
   end
   val=k{2,i};
   s=addFieldNames(s,names,val);         
end
end
function s=addFieldNames(s,names,val)
   if (numel(names)>1)
      if (~(isfield(s,names{1})))
         s.(names{1})=[];
      end
      sv=addFieldNames(s.(names{1}),{names{2:end}},val);
      s.(names{1})=sv;
   else
      if (~(isfield(s,names{1})))
         s.(names{1})=val;
      else
         if(~iscell(s.(names{1})))
            s.(names{1})={s.(names{1}),val};
         else
            s.(names{1})={s.(names{1}){:},val};
         end
      end
      
   end
end

