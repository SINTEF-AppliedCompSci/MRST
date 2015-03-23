%clear all
gravity on;
debug=false;
verbose=false;
test_hysteresis=true;
mtol=1e-2;
dmtol=3e-2;
H=10;
samples=1000;
n_test=99;
fac=3/4;
hs=linspace(0.2*H,H*(fac*0.88),n_test)';
Gt=struct('cells',struct('H',H*ones(n_test,1)));
rock=struct('perm',100*milli*darcy*ones(n_test,1),'poro',0.1*ones(n_test,1));

% make initail fluid
fluid=[];
l_fac_o=1e-9;
l_fac_g=1e-8;
fluid.surface_tension = 100;
fluid.rhoWS=1000;
fluid.rhoGS=600;
fluid.bW=@(po,varargin) (po-200*barsa)*l_fac_o+1;
fluid.BW=@(po,varargin) 1./fluid.bW(po,varargin);
fluid.bG=@(pg,varargin) (pg-200*barsa)*l_fac_g+1;
fluid.BG=@(pg,varargin) 1./fluid.bG(pg,varargin);
fluid.muW=@(po,varargin) 0.4e-3*(po*0+1);
fluid.muG=@(pg,varargin) 1e-4*(pg*0+1);
fluid_org=fluid;
l_fac_o=0e-9;
l_fac_g=0e-8;
fluid.bW=@(po,varargin) (po-200*barsa)*l_fac_o+1;
fluid.BW=@(po,varargin) 1./fluid.bW(po,varargin);
fluid.bG=@(pg,varargin) (pg-200*barsa)*l_fac_g+1;
fluid.BG=@(pg,varargin) 1./fluid.bG(pg,varargin);
fluid_incomp_org=fluid;

fluid_names = makeVEFluidsForTest();
%all_fluids= {all_fluids{}}
dsrw=linspace(0,0.2,3);
if(test_hysteresis)
    dsrg=linspace(0,0.2,3);
    %dsrw=linspace(0,0.2,3);
else
    dsrg=0;
    dsrw=0;
end
all_fluids=cell(size(fluid_names));
for j=1:numel(all_fluids)
    fluid_comp =  makeVEFluidsForTest(fluid_org, fluid_names{j},...
        'res_gas', 0,...
        'res_water', 0,...
        'Gt', Gt, 'rock', rock);
    fluid_incomp =  makeVEFluidsForTest(fluid_incomp_org, fluid_names{j},...
        'res_gas', 0,...
        'res_water', 0,...
        'Gt', Gt, 'rock', rock);
    all_fluids{j}=fluid_comp;
    all_fluids{j}.ok=true;
    ok_vector=true(10,1);
    derive_ok_vector=true(10,1);
    all_fluids{j}.ok_vector=ok_vector;
    all_fluids{j}.ok_vector=derive_ok_vector;
    for res_gas= dsrg
        for res_water= dsrw
            for kk=1:2
                if(kk==1)
                    fluid = fluid_incomp;
                else
                    fluid = fluid_comp;
                end
                fluid =  makeVEFluidsForTest(fluid, all_fluids{j}.name,...
                    'res_gas', res_gas,...
                    'res_water', res_water,...
                    'Gt', Gt, 'rock', rock);
                fluid_ok=true;
                derive_ok=true;
                for p=[233,444]*barsa;
                    s=nan(size(hs));
                    s_max=nan(size(hs));
                    kr=nan(size(hs));
                    pc=nan(size(hs));
                    for i=1:numel(hs)
                        %hs_max=hs(end)/fac;
                        hs_max=10;
                        kscale=sqrt(rock.perm(1)/rock.poro(1)).*fluid.surface_tension;
                        [ss,ppc,kkr,SH,krH, ss_max]= veRelpermTester(hs(i), p, fluid,...
                            Gt.cells.H(i), 'samples',samples,'hs_max',hs_max,'kscale',kscale);
                        s(i)=ss;
                        s_max(i)=ss_max;
                        pc(i)=ppc;
                        kr(i)=kkr;
                    end
                    ss=initVariablesADI(s);
                    akrg = fluid.krG(ss, p,'sGmax', s_max);                                      
                    dkrg=diag(akrg.jac{1});
                    krg=akrg.val;
                    dkrgnum=(krg(3:end)-krg(1:end-2))./(s(3:end)-s(1:end-2));
                    %if(any((dkrg(2:end-1)-dkrgnum)./dkrg(2:end-1)>dmtol))
                    if(any((dkrg(2:end-1)-dkrgnum)>dmtol))    
                        derive_ok=false;
                    end
                    apcog = fluid.pcWG(ss,  p,'sGmax', s_max);
                    dpcog=diag(apcog.jac{1});
                    pcog=apcog.val;
                    dpcognum=(pcog(3:end)-pcog(1:end-2))./(s(3:end)-s(1:end-2));
                    if(any((dpcog(2:end-1)-dpcognum)./dpcog(2:end-1)>dmtol))
                    %if(any((dpcog(2:end-1)-dpcognum)./(pcog(end)-pcog(1))>dmtol))                        
                        derive_ok=false;
                    end                    
                    
                    if(any(abs(krg-kr)>mtol))
                        if(verbose)
                            fprintf('Test of fluid : %s for kr failed\n', fluid.name);
                            %fprintf('veRelpermTester say : %d\n', kr);
                            %fprintf('Calculated value is say : %d\n', krg);
                        end
                        fluid_ok=false;
                        
                    end
                    if(any(abs(pcog-pc)>mtol*pc(end)))
                        if(verbose)
                            fprintf('Test of fluid : %s for pc failed\n', fluid.name);
                            fprintf('Pressure is %d\n',p/barsa)
                            fprintf('%d %d\n',[pcog';pc'])
                            %fprintf('veRelpermTester say : %d\n', kr);
                            %fprintf('Calculated value is say : %d\n', krg);
                        end
                        fluid_ok=false;
                    end
                    if(~(fluid_ok && derive_ok))
                        if(debug)
                            %%
                            figure(1),clf
                            disp(all_fluids{j}.name)
                            subplot(1,2,1)
                            plot(s,krg,'*',s,kr,[s;1],[s;1])
                            subplot(1,2,2)
                            plot(s,pcog,'*',s,pc)
                            return
                        end
                    end
                end
                if(~fluid_ok)
                    ok_vector(1)=false;
                    if(kk==1)
                        ok_vector(2) = false;
                        if((res_gas==0) && (res_water==0))
                            ok_vector(3) = false;
                        end
                    else
                        ok_vector(4)=false;
                        if((res_gas==0) && (res_water==0))
                            ok_vector(5) = false;
                        end
                    end
                end
                if(~derive_ok)
                    derive_ok_vector(1)=false;
                    if(kk==1)
                        derive_ok_vector(2) = false;
                        if((res_gas==0) && (res_water==0))
                            derive_ok_vector(3) = false;
                        end
                    else
                        derive_ok_vector(4)=false;
                        if((res_gas==0) && (res_water==0))
                            derive_ok_vector(5) = false;
                        end
                    end
                end
            end
        end
    end
    all_fluids{j}.ok=fluid_ok;
    all_fluids{j}.ok_vector=ok_vector;
    all_fluids{j}.derive_ok_vector=derive_ok_vector;
end
fprintf('Test of fluid : %20s', 'name')
%for k=1:5
%    fprintf('\t %5s',)
%end
fprintf('\t %15s','all')
fprintf('\t %15s','incomp')
fprintf('\t %15s','incomp_nores')
fprintf('\t %15s','comp')
fprintf('\t %15s','comp_nores')
fprintf('\n')
for j=1:numel(all_fluids)
    fprintf('Test of fluid values: %20s', all_fluids{j}.name)
    for k=1:5
    if(all_fluids{j}.ok_vector(k))
        fprintf('\t %15s','OK');
    else
        fprintf('\t %15s','Failed');
    end
    end
    fprintf('\n')
end
fprintf('\n')
fprintf('Test of fluid derivative: %20s', 'name')
%for k=1:5
%    fprintf('\t %5s',)
%end
fprintf('\t %15s','all')
fprintf('\t %15s','incomp')
fprintf('\t %15s','incomp_nores')
fprintf('\t %15s','comp')
fprintf('\t %15s','comp_nores')
fprintf('\n')
for j=1:numel(all_fluids)
    fprintf('Test of fluid : %20s', all_fluids{j}.name)
    for k=1:5
    if(all_fluids{j}.derive_ok_vector(k))
        fprintf('\t %15s','OK');
    else
        fprintf('\t %15s','Failed');
    end
    end
    fprintf('\n')
end
