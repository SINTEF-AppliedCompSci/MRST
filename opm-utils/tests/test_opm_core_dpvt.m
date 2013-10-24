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
my_tol=1e+2;
pmin=200.1;pmax=300;npress=1000;
param.output_dir='data/test_output';
param.tab_size_kr=3;
param.dead_tab_size=40;
param.input_filename='data/test_pressure.txt';
param.relperm_output='data/relperm.txt';
param.matrix_output='data/matrix.txt';
param.mu_output='data/c.txt';
param.rs_output='data/rs.txt';
param.b_output='data/b.txt';
%layers=[1:85];
have_gas=true;   % test with gas
nw=3;alpha_w=3;% points,exponent SWOF
ng=3;alpha_g=2;% points,exponent GSOF
clear deck;
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
deck.PROPS.PVTW=[200 1 0.0 0.5 0];

%oil_pvt='PVCDO';gas_pvt='PVDG';
oil_pvt='PVDO';gas_pvt='PVDG';

% always gas gas and oil
%deck.SOLUTION.SOIL(:)=1-(initsat+initgas);
%deck.SOLUTION.SWAT(:)=initsat;
switch oil_pvt
   case 'PVCDO'      
      deck.PROPS.PVCDO=[200 1 0 5 0];
      %deck.PROPS.PVTW(:,3)=0.5e-5*barsa/psia;
      deck.PROPS.PVCDO(:,3)=0.5e-3*barsa/psia;
      is_deadoil=true;
   case 'PVDO'       
      %deck.PROPS.PVDO{1}(:,2)=deck.PROPS.PVDO{1}(1,2);
      %deck.PROPS.PVTW(1,4)=deck.PROPS.PVDO{1}(1,end)
      %deck.PROPS.PVDO{1}(:,end)=deck.PROPS.PVTW(1,4);
      deck.PROPS.PVDO{1}=[50 0.9  2;
                       100 0.7  2;
                       600 0.1  2];
      deck.PROPS.PVDO{1}=[100 0.7  2;
                       600 0.1  2];              
      is_deadoil=true;
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
                 600 0.4  10];
   otherwise
      error('No case')
   end
end
writeDeck(deck,fullfile('data','test_deck'))
param.deck_filename=fullfile('data','test_deck','test_deck.DATA');
% define saturations to use
names = {'water', 'oil', 'gas'};              
phase = isfield(deck.RUNSPEC, upper(names));
np=sum(phase);
fid=fopen(param.input_filename,'w');
press=linspace(pmin,pmax,npress)*barsa;
is_dead_oil=true;
if(is_dead_oil)
   for i=1:numel(press) 
     fprintf(fid,'%f ',press(i));
     for j=1:np-1
        fprintf(fid,'1.0 ');
     end
     fprintf(fid,'1.0\n');
   end
else
   error();
end
fclose(fid);
%% write param file
!input_filename=fullfile('data','test_run.param');
paramStructToParamFile(param,input_filename);
%{
executable='pvt_test';
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH;';
mydir=opm_code_dir();
if(false)
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

executable='tests/pvt_test';
mydir=opm_dir(myself);
command = [simcommand(mydir, executable), ' ', input_filename]

a=system(command);
if(a>0)
   error(['Command failed : ', command]);
end
%% define fluid in matlab
deck=readEclipseDeck(param.deck_filename);
deck=convertDeckUnits(deck);
fluid=initEclipseFluid(deck);
% defind 
res.opm=[];
res.opm.B=load(param.b_output);
res.opm.dB=res.opm.B(:,np+1:end);
res.opm.B=res.opm.B(:,1:np);
res.opm.RS=load(param.rs_output);
res.opm.dRS=res.opm.RS(:,np+1:end);
res.opm.RS=res.opm.RS(:,1:np);
Ao=load(param.matrix_output);
res.opm.A=Ao(:,1:np*np);
res.opm.dA=Ao(:,np*np+1:end);
res.opm.invB=1./res.opm.B;
res.opm.dinvB=-res.opm.dB./res.opm.B.^2;
res.opm.c=-res.opm.dB./res.opm.B;
pz=load(param.input_filename);
p=pz(:,1);
z=pz(:,2:1+np);

[kr.mrst,dkr.mrst]=fluid.pvt(p,z);
[res.mrst.c, res.mrst.rho,.... 
   res.msrt.mu, res.mrst.u, res.mrst.RS,...
   res.mrst.B,res.mrst.A, res.mrst.dA, res.mrst.dB, res.mrst.dRS] = fluid.pvt(p, z);
res.mrst.B=reshape(full(diag(res.mrst.B)),np,[])';
res.mrst.dB=reshape(full(diag(res.mrst.dB)),np,[])';
res.mrst.invB=1./res.mrst.B;
res.mrst.dinvB=-res.mrst.dB./res.mrst.B.^2;
diff_B=res.mrst.B-res.opm.B;
var={'RS','A','dRS','dA'};
for d=1:numel(var)
   tmp=zeros(numel(p),np*np);
   for i=1:numel(p)
      a=full(res.mrst.(var{d})((i-1)*np+1:i*np,(i-1)*np+1:i*np));
      tmp(i,:)=a(:);
   end
   res.mrst.(var{d})=tmp;
end
figure(33),clf
if(is_dead_oil)
   nf=4;
   subplot(nf,1,1),cla
   plot(p(:,1)/barsa,[res.mrst.B,res.opm.B],'*-')
   legend('msrt 1','msrt 2','mrst 3','opm 1','opm 2','opm 3')
   title('B factor')
   subplot(nf,1,2),cla
   plot(p(:,1)/barsa,[res.mrst.dB,res.opm.dB],'*-')
   legend('msrt 1','msrt 2','mrst 3','opm 1','opm 2','opm 3')
   title('dB')
   subplot(nf,1,3),cla   
   plot(p(:,1)/barsa,[res.mrst.dinvB,res.opm.dinvB],'*-')
   legend('msrt 1','msrt 2','mrst 3','opm 1','opm 2','opm 3')
   title('dinvB')
   subplot(nf,1,4),cla   
   plot(p(:,1)/barsa,[res.mrst.invB,res.opm.invB],'*-')
   legend('msrt 1','msrt 2','mrst 3','opm 1','opm 2','opm 3')
   title('invB')
else

end


%diff_dkr=dkr.mrst-dkr.opm;
%disp('Difference in B')
%disp(diff_B)
%% Difference between opm and mrst is since opm resample 1/B and then
%% assume d(1/B)dp=C in the new table which is the same as the above only
%% if all points in the initial table is in the last table.
if(is_dead_oil)
   mycases={'mrst','opm'};
   var={'B','RS','A','invB'}
   dvar={'dB','dRS','dA','dinvB'}
   for d=1:numel(var)
      for k=1:numel(mycases)
         max_error=zeros(numel(p),size(res.(mycases{k}).(var{d})(1,:),2));
         for i=1:numel(p)-1;
            dp=p(i+1)-p(i);
            dB_num=res.(mycases{k}).(var{d})(i+1,:)-res.(mycases{k}).(var{d})(i,:);
            dB_ana=res.(mycases{k}).(dvar{d})(i,:).*dp;
            rel_err=abs((dB_num-dB_ana)/max(abs(dB_ana)));
            if(any(rel_err>my_tol) && ~all((dB_num-dB_ana)==0))
               error([var{d}, ' failed for ',mycases{k}])
            end
            max_error(i,:)=rel_err;
         end
         disp([var{d} ,' success for ',mycases{k} ,' max relative error ', num2str(max(max_error))])
      end
   end          
else
 error('Live oil check not implemented')   
end


%%

