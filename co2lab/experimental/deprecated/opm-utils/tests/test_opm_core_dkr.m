mrstModule add opm-utils
mrstModule('add',fullfile(ROOTDIR,'mex','libgeometry'))
mrstModule add spe10
mrstModule add eclipse
clear param;
param.output_dir='data/test_output';
param.input_filename='data/test_saturations.txt';
param.relperm_output='data/relperm.txt';
%param.sat_tab_size='0';
%%
dkr_tol=1e-2;
have_gas=true;   % test with gas
nw=201;alpha_w=3;% points,exponent SWOF
ng=201;alpha_g=2;% points,exponent GSOF
sw=linspace(0,1,100); %test grid sw
if(have_gas)
   sg=linspace(0,1,100);% test grid in sg
else
   sg=0;
end
%% set up deck
deck=deckQFS([1 1 1],[1 1 1]);
nstep=1;
%deck.PROPS.STONE=1;
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

%
S=zeros(numel(sw)*numel(sg),sum(phase));
k=0;
fid=fopen(param.input_filename,'w');
for i=1:numel(sw)
   if(have_gas)
      for j=1:numel(sg);        
         if((sw(i)+sg(j))<1)
           k=k+1;
           I(k)=i;J(k)=j;SW(k)=sw(i);SG(k)=sg(j);SO(k)=1-(sw(i)+sg(j));
           S(k,:)=[sw(i),1-(sw(i)+sg(j)),sg(j)];
           K(k)=k;
           fprintf(fid,'%f ',S(k,:));
           fprintf(fid,'\n');
         end
      end
   else
      k=k+1;
      I(k)=i;J(k)=j;SW(k)=sw(i);So(k)=1-sw(i);K(k)=k;
      S(k,:)=[1-sw(i),sw(i)];
      fprintf(fid,'%f ',S(k,:));
      fprintf(fid,'\n');
   end 
end
fclose(fid)
S=S(1:k,:);
%save(param.input_file,'--ASCII',S);
%% write param file
input_filename=fullfile('data','test_run.param');
paramStructToParamFile(param,input_filename);

%{
executable='relperm_test';
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH;'
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
%}

% New style command creation.
executable='tests/relperm_test';
mydir=opm_dir(myself);
command = [simcommand(mydir, executable), ' ', input_filename]

a=system(command)
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

%% Compare forward differences with analytic estimates
dsw = 0.0;
dsg = 0.0;
if(np==3)
  S_W=sparse(I,J,SW);
  S_O=sparse(I,J,SO);
  S_G=sparse(I,J,SG);
  K_K=sparse(I,J,K);
  mytypes={'mrst','opm'};
  %mytypes={'mrst'}
  for kk=1:numel(mytypes)
     for i=1:numel(sw)-1;
        for j=1:numel(sg)-1;
           if((sw(i+1)+sg(j+1))<(1-dsw-dsg))
              k_0=K_K(i,j);Soo=[S_W(i,j),S_O(i,j),S_G(i,j)];
              % x direction
              k_1=K_K(i+1,j);Sxo=[S_W(i+1,j),S_O(i+1,j),S_G(i+1,j)];
              dS=Sxo-Soo;
              dkr_num=kr.(mytypes{kk})(k_1,:)-kr.(mytypes{kk})(k_0,:);
              dkr_ana=reshape(dkr.(mytypes{kk})(k_0,:),np,np)*dS';
              if(norm(dkr_num'-dkr_ana,'inf')>dkr_tol)
                 error(['Error S1 direction ',mytypes{kk}])
              end
              k_2=K_K(i,j+1);Soy=[S_W(i,j+1),S_O(i,j+1),S_G(i,j+1)];
              dS=Soy-Soo;
              dkr_num=kr.(mytypes{kk})(k_2,:)-kr.(mytypes{kk})(k_0,:);
              dkr_ana=reshape(dkr.(mytypes{kk})(k_0,:),np,np)*dS';
              if(norm(dkr_num'-dkr_ana,'inf')>dkr_tol)
                 error(['Error S2 direction ',mytypes{kk}])
              end
              k_3=K_K(i+1,j+1);Syy=[S_W(i+1,j+1),S_O(i+1,j+1),S_G(i+1,j+1)];
              dS=Syy-Soo;
              dkr_num=kr.(mytypes{kk})(k_3,:)-kr.(mytypes{kk})(k_0,:);
              dkr_ana=reshape(dkr.(mytypes{kk})(k_0,:),np,np)*dS';
              if(norm(dkr_num'-dkr_ana,'inf')>dkr_tol)
                 error(['Error diagonal ',mytypes{kk}])
              end
           end
        end
     end
     disp(['Success for ', mytypes{kk}])
  end
elseif(np==2)
  for i=1:size(S,1)-1;
    dS=S(i+1,:)-S(i,:);
    %dS=[dS,-dS];    
    % test matlab
    mytypes={'mrst','opm'};
    for k=1:numel(mytypes)
       dk_num=kr.(mytypes{k})(i+1,:)-kr.(mytypes{k})(i,:);
       dk_ana=reshape(dkr.(mytypes{k})(i,:),2,2)*dS';
       sprintf('%f ',dk_num-dk_ana')';
       sprintf('\n');
       if(norm(dk_num'-dk_ana,'inf')>dkr_tol)
          error('Twophase test failed');
       end
    end
  end
end

%% Compare MRST and OPM kr to each other
krdiff = max(max(kr.opm - kr.mrst))
if krdiff < 1e-4
    disp('OPM and MRST kr are close');
else
    error('OPM and MRST kr differ.');
end
