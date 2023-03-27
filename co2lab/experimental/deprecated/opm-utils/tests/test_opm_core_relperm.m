mrstModule add opm-utils
mrstModule('add',fullfile(ROOTDIR,'mex','libgeometry'))
mrstModule add spe10
clear param;

if ~isdir('data'),
   [succ, msg, id] = mkdir('data');

   if ~succ,
      error(id, 'Failed to create directory data: %s', msg);
   end
end

param.output_dir='data/test_output';
param.tab_size=3;
param.input_filename='data/test_saturations.txt';
param.relperm_output='data/relperm.txt';
%layers=[1:85];
have_gas=true;   % test with gas
nw=3;alpha_w=3;% points,exponent SWOF
ng=3;alpha_g=2;% points,exponent GSOF
deck=deckQFS([1 1 1],[1 1 1]);
nstep=1;
%deck.PROPS.STONE1=1;
deck.PROPS.SIMPLE=1;
deck.SCHEDULE.step.control=ones(nstep,1);
deck.SCHEDULE.step.val=ones(nstep,1)*2000/nstep;
%% define fluid
myfields={'PVCDO','PVTW'};
for i=1:numel(myfields)
   deck.PROPS=rmfield(deck.PROPS,myfields{i});
end
deck.PROPS.PVTW=[0 1 0 0.5 0];

%oil_pvt='PVCDO';gas_pvt='PVDG';
oil_pvt='PVDO';gas_pvt='PVDG';

% always gas gas and oil
%deck.SOLUTION.SOIL(:)=1-(initsat+initgas);
%deck.SOLUTION.SWAT(:)=initsat;
switch oil_pvt
   case 'PVCDO'      
      deck.PROPS.PVCDO=[0 1 0 5 0];
      %deck.PROPS.PVTW(:,3)=0.5e-5*barsa/psia;
      deck.PROPS.PVCDO(:,3)=0.5e-5*barsa/psia;
   case 'PVDO'       
      %deck.PROPS.PVDO{1}(:,2)=deck.PROPS.PVDO{1}(1,2);
      %deck.PROPS.PVTW(1,4)=deck.PROPS.PVDO{1}(1,end)
      %deck.PROPS.PVDO{1}(:,end)=deck.PROPS.PVTW(1,4);
      deck.PROPS.PVDO{1}=[0 1  2;
                       100 1  2;
                       600 1  2;
                       1000 0.4  2];
   otherwise
      error('No case')
end
s=linspace(0,1,nw)';
deck.PROPS.SWOF{1}=[s,s.^alpha_w,(1-s).^alpha_w,s*0];
if(have_gas)
   initgas=0.001;
   s=linspace(0,1,ng)';
   deck.PROPS.SGOF{1}=[s,s.^alpha_g,(1-s).^alpha_g,s*0];
   deck.RUNSPEC.GAS=1;
   deck.SOLUTION.SGAS=ones(size(deck.SOLUTION.SWAT))*initgas;
   switch gas_pvt
      case 'PVTG'
         error('not implemented')
      case 'PVDG'
      %deck.PROPS.PVDO{1}(:,end)=deck.PROPS.PVTW(1,4);
      deck.PROPS.PVDG{1}=[100 1  10;
                 600 0.6  10];
   otherwise
      error('No case')
   end
end
writeDeck(deck,fullfile('data','test_deck'))
param.deck_filename=fullfile('data','test_deck','test_deck.DATA');

% define saturations to use
names = {'water', 'oil', 'gas'};              
phase = isfield(deck.RUNSPEC, upper(names));
sw=linspace(0,1,100);
if(have_gas)
   sg=linspace(0,1,100);
else
   sg=0;
end
S=zeros(numel(sw)*numel(sg),sum(phase));
k=0;
fid=fopen(param.input_filename,'w');
for i=1:numel(sw)
   if(have_gas)
      for j=1:numel(sg);        
         if((sw(i)+sg(j))<=1)
            k=k+1;            
            S(k,:)=[sw(i),1-(sw(i)+sg(j)),sg(j)];
            fprintf(fid,'%f ',S(k,:));
            fprintf(fid,'\n');
         end
      end
   else
      k=k+1;
      S(k,:)=[1-sw(i),sw(i)];
      fprintf(fid,'%f ',S(k,:));
      fprintf(fid,'\n');
   end 
   
end
S=S(1:k,:);
fclose(fid)
%save(param.input_file,'--ASCII',S);
%% write param file
input_filename=fullfile('data','test_run.param');
paramStructToParamFile(param,input_filename);
executable='relperm_test';
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH;';
mydir=opm_code_dir();
if(true)
   if(~isempty(param))
      command=[ldfix, fullfile(mydir,'builds/db/opm-core/tests/',executable),' ',input_filename];
   else
      command=[ldfix, fullfile(mydir,'builds/db/opm-core/tests/',executable)];
   end
else
   if(~isempty(param))
      command=[ldfix, fullfile(mydir,'qtcreator-build/',executable),' ',input_filename];
   else
      command=[ldfix, fullfile(mydir,'qtcreator-build/',executable)];
   end
end
a=system(command);
if(a>0)
   disp(['Command failed : ', command]);
end
%% define fluid in matlab
fluid=initEclipseFluid(deck);
skr=load(param.relperm_output);
np=sum(phase);
s=skr(:,1:np);
kr.opm=skr(:,np+1:2*np);
dkr.opm=skr(:,2*np+1:2*np+np*np);
[kr.mrst,dkr.mrst]=fluid.relperm(s);

diff_kr=kr.mrst-kr.opm;
diff_dkr=dkr.mrst-dkr.opm;
disp('Difference in kr')
disp(diff_kr)
%disp('dkr diff')
%disp(diff_dkr)
%%
figure(33),clf
if(np==2)
   plot(s(:,1),[kr.mrst,kr.opm])
   legend('msrt 1','msrt 2','opm 1','opm 2')
else
   for i=1:np
      subplot(np,1,i),cla
      plot(s(:,i),[kr.mrst,kr.opm],'*')
   end
   legend('msrt 1','msrt 2','mrst 3','opm 1','opm 2','opm 3')
end

%%
myfields={'mrst','opm'};
if(np==2) 
  % use dS1 as variable
  for i=1:numel(myfields)
     myfield=myfields{i};
     krdS1.(myfield)=[dkr.(myfield)(:,1)- dkr.(myfield)(:,3),dkr.(myfield)(:,2)- dkr.(myfield)(:,4)];
  end
 
 disp('Difference in dkr restrikted')
 disp(krdS1.mrst-krdS1.opm)
elseif(np==3)
  % use dS1 and dS2
  for i=1:numel(myfields)
     myfield=myfields{i};
     krdS1.(myfield)=[dkr.(myfield)(:,1)- dkr.(myfield)(:,7) , dkr.(myfield)(:,2)- dkr.(myfield)(:,8),dkr.(myfield)(:,3)- dkr.(myfield)(:,9)];
     krdS2.(myfield)=[dkr.(myfield)(:,4)- dkr.(myfield)(:,7),dkr.(myfield)(:,5)- dkr.(myfield)(:,8),dkr.(myfield)(:,6)- dkr.(myfield)(:,9)];
  end
  disp('Difference in dkrdS1 restrikted')
  disp(krdS1.mrst-krdS1.opm)
  disp('Difference in dkrdS1 restrikted')
  disp(krdS2.mrst-krdS2.opm)
else
  error('Only two and three phase is valid');
end
%% check derivatives numerically;
figure(44)
if(np==3)
   [S1,S2] = MESHGRID(linspace(0,1,100),linspace(0,1,100));
   for i=0:np-1;
      Z_mrst=griddata(s(:,1),s(:,2),kr.mrst(:,i+1),S1,S2,'linear');
      Z_mrst(S1(:)+S2(:)>1)=nan;
      subplot(np,2,2*i+1)
      mesh(S1,S2,Z_mrst);
      title(['mrst kr',num2str(i)]);
      xlabel('S1')
      ylabel('S2')
      view([0.1 -1 1])
      Z_opm=griddata(s(:,1),s(:,2),kr.opm(:,i+1),S1,S2,'linear');
      Z_opm(S1(:)+S2(:)>1)=nan;
      subplot(np,2,2*i+2)
      mesh(S1,S2,Z_opm);
      title(['opm kr',num2str(i)])
      xlabel('S1')
      ylabel('S2')
      view([0 1 1])
      view([0.1 -1 1])
      for k=1:2;
         subplot(np,2,2*i+k)
         if(i==0)
           view([-0.3 -1 1]) 
         elseif(i==1)
            view([1 -0.1 1])
         elseif(i==2)
            view([0.1 1 1])
         end
      end
   end
end
   
