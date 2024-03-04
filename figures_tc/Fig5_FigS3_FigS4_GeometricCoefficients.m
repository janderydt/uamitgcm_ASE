function Fig5_FigS3_FigS4_GeometricCoefficients

addpath(getenv("froot_tools"));

endtime = datenum('01012215','ddmmyyyy');

runID = {'PTDC_002','PTDC_003'};
runIDlinestyle = {'-','-','-','-'};
runIDlinewidth = [2 2];
basins = {'PIG','TW','CR','DT'};
basinslegend = {'Pine Island','Thwaites','Crosson','Dotson'};
%color = [[51 34 136];[17 119 51];[221 204 119];[136 34 85]]/255;
%color = [[51 34 136];[17 119 51];[133 87 35];[136 34 85]]/255;
%color=[[166,206,227];[31,120,180];[178,223,138];[51,160,44]]/255;
%color = [[217,95,2];[117,112,179];[231,41,138];[102,166,30]]/255;
%color = [[228,26,28];[55,126,184];[77,175,74];[152,78,163]]/255;
%colorsoft = [[251,180,174];[179,205,227];[204,235,197];[222,203,228]]/255;
color = [[31,120,180];[27,158,119];[246,190,0];[152,78,163]]/255;
colorsoft = [[166,206,227];[141,211,199];[251,219,101];[190,174,212]]/255;

rhoConst = 1024;
Cp = 3974;
Lf = 334e3; %J/kg
Cd = 0.0025;
GammT = 0.015;
m0 = 2.2e-4;

ylabels = {'$\widetilde{m}$','$\widetilde{T}_{\star {\rm IF}}{}^2$',...
    '$\widetilde{\mu}^2$','$\widetilde{\epsilon}_T$','$\widetilde{\epsilon}_U$'};
panellabel = {'a','b','c','d','e','f','g','h','i','j'};

H1=fig('units','inches','width',52.5*12/72.27,'height',60*12/72.27,'fontsize',10,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);

tlo_fig = tiledlayout(5,numel(runID),'TileSpacing','tight');
for i = 1:5*numel(runID)
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

H2 = fig('units','inches','width',40*12/72.27,'height',40*12/72.27,'fontsize',12,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);

tlo_fig2 = tiledlayout(2,2,'TileSpacing','tight');
for i = 1:4
    ax_fig2(i) = nexttile(tlo_fig2); hold on;
end

H3 = fig('units','inches','width',40*12/72.27,'height',40*12/72.27,'fontsize',12,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);

tlo_fig3 = tiledlayout(2,2,'TileSpacing','tight');
for i = 1:4
    ax_fig3(i) = nexttile(tlo_fig3); hold on;
end
% 
% H3 = fig('units','inches','width',100*12/72.27,'height',55*12/72.27,'fontsize',16,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);
% 
% tlo_fig3 = tiledlayout(3,3,'TileSpacing','tight');
% ax_fig3(1) = nexttile(tlo_fig3,[3 1]); hold on;
% for i = 2:5
%     ax_fig3(i) = nexttile(tlo_fig3,[1 1]); hold on;
% end
% ax_fig3(6) = nexttile(tlo_fig3,[1 2]); hold on;
% 
% H4 = fig('units','inches','width',80*12/72.27,'height',35*12/72.27,'fontsize',16,'font','Helvetica');
% 
% tlo_fig4 = tiledlayout(1,3,'TileSpacing','tight');
% for i = 1:3
%    ax_fig4(i) = nexttile(tlo_fig4); hold on;
% endf

for jj=1:numel(runID)

    load(['HeatVolumeTransport_IceFront_below400m_PTDC_001.mat']);
    MITTime_IF_001 = MITTime;
    section_IF_001 = section;

    load(['HeatVolumeTransport_moving400mdraft_PTDC_001.mat']);
    MITTime_001 = MITTime;
    section_001 = section;
    integral2D_001 = integral2D;

    load(['HeatVolumeTransport_IceFront_below400m_',runID{jj},'.mat']);
    MITTime_IF = MITTime;
    section_IF = section;

    load(['HeatVolumeTransport_moving400mdraft_',runID{jj},'.mat']);
       
    SF = load('CavityStreamFunctions_PTDC_002_monthly.mat');

    for gg=[numel(basins):-1:1]

        basin = basins{gg};
        switch basin
            case 'PIG'
                xpoly = [-1699e3 -1699e3 -1530e3 -1530e3];
                ypoly = [-370e3 -220e3 -220e3 -370e3];
            case 'TW'
                xpoly = [-1620e3 -1620e3 -1500e3 -1500e3];
                ypoly = [-520e3 -370e3 -370e3 -520e3];
            case 'CR'
                xpoly = [-1630 -1485 -1450 -1450 -1630]*1e3;
                ypoly = [-600 -657 -657 -520 -520]*1e3; 
            case 'DT'
                xpoly = [-1630 -1485 -1450 -1450 -1630]*1e3;
                ypoly = [-600 -657 -657 -700 -700]*1e3; 
            case 'AS'
                xpoly = [-1700e3 -1700e3 -1450e3 -1450e3];
                ypoly = [-700e3 -220e3 -220e3 -700e3];
            case ''
                xpoly = [-1700e3 -1700e3 -1450e3 -1450e3];
                ypoly = [-700e3 -200e3 -200e3 -700e3];
        end
        
        %% spinup
        [~,end_001] = min(abs(MITTime_001 - datenum('01122015','ddmmyyyy')));
        [~,end_IF_001] = min(abs(MITTime_IF_001 - datenum('01122015','ddmmyyyy')));

        TminTf_IF=[]; kk=1;
       
        for tt=1:end_IF_001

            I = find(inpoly([section_IF_001(tt).xmid(:) section_IF_001(tt).ymid(:)],[xpoly(:) ypoly(:)]));
            
            % thermal driving, heat flux, volume flux, cross section
            TminTf_IF_001(kk) = mean(section_IF_001(tt).monthly.TminTf_belowz_in_depthmean(I),'all','omitnan');
            %TminTf_IF_001(kk) = mean(section_IF_001(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');
            kk=kk+1;
            
       end
       
       MITTime_IF_001 = MITTime_IF_001(1:end_IF_001);
       [MITTime_IF_001,ia,ic] = unique(MITTime_IF_001);
       TminTf_IF_001 = TminTf_IF_001(ia);

       TminTf_001=[]; Ustar_bl_001=[]; TminTf_bl_001=[]; melt_001=[]; kk=1;

       for tt=1:end_001

            I = find(inpoly([section_001(tt).xmid(:) section_001(tt).ymid(:)],[xpoly(:) ypoly(:)]));
            
            % thermal driving, heat flux, volume flux, cross section
            TminTf_001(kk) = mean(section_001(tt).monthly.TminTf_belowz_in_depthmean(I),'all','omitnan');
            %TminTf_001(kk) = mean(section_001(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');
            Ustar_bl_001(kk) = integral2D_001.(basin).monthly.blUStar_mean(tt);
            TminTf_bl_001(kk) = integral2D_001.(basin).monthly.TminTf_mean(tt);
            melt_001(kk) = integral2D_001.(basin).monthly.melt_integral(tt)./integral2D_001.(basin).monthly.ISarea_integral(tt)*1e12/1e3; %m/yr

            kk= kk+1;
    
       end

       MITTime_001 = MITTime_001(1:end_001);
       [MITTime_001,ia,ic] = unique(MITTime_001);
       TminTf_001 = TminTf_001(ia); 
       Ustar_bl_001 = Ustar_bl_001(ia);
       TminTf_bl_001 = TminTf_bl_001(ia);
       melt_001 = melt_001(ia);

       mu_001 = TminTf_001./interp1(MITTime_IF_001,TminTf_IF_001,MITTime_001);
       epsU_001 = Ustar_bl_001./TminTf_001;
       epsT_001 = TminTf_bl_001./TminTf_001;


        %% transient
       if contains(runID{jj},{'PTDC_001','PTDC_000'})
            [~,start_MITTime] = min(abs(MITTime - datenum('01011995','ddmmyyyy')));
            [~,start_MITTime_IF] = min(abs(MITTime_IF - datenum('01011995','ddmmyyyy')));
       elseif contains(runID{jj},{'PTDC_002','PTDC_003','PTDC_004'})
            [~,start_MITTime] = min(abs(MITTime - datenum('01012015','ddmmyyyy')));
            [~,start_MITTime_IF] = min(abs(MITTime_IF - datenum('01012015','ddmmyyyy')));
       end
       
       [~,end_MITTime] = min(abs(MITTime - endtime));
       [~,end_MITTime_IF] = min(abs(MITTime_IF - endtime));

       TminTf_IF=[]; overturning_IF=[]; kk=1;
       TminTf_IF_depthmean=[]; TminTf_IF_belowz_depthmean=[]; TminTf_IF_in_depthmean=[];
       
       for tt=start_MITTime_IF:end_MITTime_IF

            I = find(inpoly([section_IF(tt).xmid(:) section_IF(tt).ymid(:)],[xpoly(:) ypoly(:)]));
            
            % thermal driving, heat flux, volume flux, cross section
            TminTf_IF_depthmean(kk) = mean(section_IF(tt).monthly.TminTf_depthmean(I),'all','omitnan');
            TminTf_IF_belowz_depthmean(kk) = mean(section_IF(tt).monthly.TminTf_belowz_depthmean(I),'all','omitnan');
            TminTf_IF_in_depthmean(kk) = mean(section_IF(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');
            TminTf_IF(kk) = mean(section_IF(tt).monthly.TminTf_belowz_in_depthmean(I),'all','omitnan');
            %overturning_IF(kk)=section_IF(tt).monthly.maxoverturning.(basin);

            kk=kk+1;
            
       end

       MITTime_IF = MITTime_IF(start_MITTime_IF:end_MITTime_IF);

       %figure(222+gg); hold on;
       %plot(MITTime_IF,TminTf_IF_depthmean,MITTime_IF,TminTf_IF_belowz_depthmean,MITTime_IF,TminTf_IF_in_depthmean,MITTime_IF,TminTf_IF);
        

       [MITTime_IF,ia,ic] = unique(MITTime_IF);
       TminTf_IF = TminTf_IF(ia);
       %overturning_IF = overturning_IF(ia);

       Ustar_bl =[]; TminTf_bl=[]; TminTf_DI=[]; melt=[]; 
       bgradx=[]; bgrady=[]; uvel=[]; vvel=[]; 
       gradxrho=[]; gradyrho=[]; uvel_bt=[]; vvel_bt=[];
       kk=1;

       for tt=start_MITTime:end_MITTime

            I = find(inpoly([section(tt).xmid(:) section(tt).ymid(:)],[xpoly(:) ypoly(:)]));
            
            % thermal driving, heat flux, volume flux, cross section
            Tin(kk) = mean(section(tt).monthly.T_belowz_in_depthmean(I),'all','omitnan');
            TminTf_DI(kk) = mean(section(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');
            %vel_in(kk) = mean(section(tt).monthly.velperp_in_depthmean(I),'all','omitnan');

            Ustar_bl(kk) = integral2D.(basin).monthly.blUStar_mean(tt);
            TminTf_bl(kk) = integral2D.(basin).monthly.TminTf_mean(tt);

            melt(kk) = integral2D.(basin).monthly.melt_integral(tt)./integral2D.(basin).monthly.ISarea_integral(tt)*1e12/1e3; %m/yr

            draft(kk) = integral2D.(basin).monthly.draft_mean(tt);
            bgrady(kk)=integral2D.(basin).monthly.bgrady_mean(tt);
            bgradx(kk)=integral2D.(basin).monthly.bgradx_mean(tt);

            vvel_bl(kk)=integral2D.(basin).monthly.blVVEL_mean(tt);
            uvel_bl(kk)=integral2D.(basin).monthly.blUVEL_mean(tt);
            vel_bl(kk) = integral2D.(basin).monthly.blVEL_mean(tt);
            
            uvel_bt(kk) = integral2D.(basin).monthly.btUVEL_mean(tt);
            vvel_bt(kk) = integral2D.(basin).monthly.btVVEL_mean(tt);
            vel_bt(kk) =  integral2D.(basin).monthly.btVEL_mean(tt);
            %gradxrho(kk)=integral2D.(basin).monthly.gradxRho_mean(tt);
            %gradyrho(kk)=integral2D.(basin).monthly.gradyRho_mean(tt);

            uvelratio(kk) = integral2D.(basin).monthly.UVELratio(tt);
            vvelratio(kk) = integral2D.(basin).monthly.VVELratio(tt);
            gradbx_abs(kk) = integral2D.(basin).monthly.gradbx_abs(tt);
            gradby_abs(kk) = integral2D.(basin).monthly.gradby_abs(tt);
            velratio_gradb(kk) = integral2D.(basin).monthly.VELratio_gradb(tt);
            epsT_new(kk) =  integral2D.(basin).monthly.epsT(tt);

            kk= kk+1;
    
       end

       MITTime = MITTime(start_MITTime:end_MITTime);
       [MITTime,ia,ic] = unique(MITTime);
       TminTf_DI = TminTf_DI(ia); 
       melt = melt(ia);
       Ustar_bl = Ustar_bl(ia);
       TminTf_bl = TminTf_bl(ia);
%        draft = draft(ia);
%        Tin = Tin(ia);
        bgradx = bgradx(ia);
        bgrady = bgrady(ia);
%        uvel_bl = uvel_bl(ia);
%        vvel_bl = vvel_bl(ia);
%        uvel_bt = uvel_bt(ia);
%        vvel_bt = vvel_bt(ia);
%        vel_bl = vel_bl(ia);
%        vel_bt = vel_bt(ia);
%        uvel_shear = uvel_shear(ia);
%        vvel_shear = vvel_shear(ia);
%        vel_shear = vel_shear(ia);
%        %vel_in = vel_in(ia);
        %gradxrho = gradxrho(ia);
        %gradyrho = gradyrho(ia);
        uvelratio = uvelratio(ia);
        vvelratio = vvelratio(ia);
        gradbx_abs = gradbx_abs(ia);
        gradby_abs = gradby_abs(ia);
        velratio_gradb = velratio_gradb(ia);
        epsT_new = epsT_new(ia);
        % 
        mu = TminTf_DI./interp1(MITTime_IF,TminTf_IF,MITTime);
        epsU = Ustar_bl./TminTf_DI;
        epsT = TminTf_bl./TminTf_DI;
% %        bsf = interp1(SF.time,SF.bsf.(basin),MITTime);
% %        osf = interp1(SF.time,SF.osf.(basin),MITTime)
% 
%        vel_bl_tmp = sqrt(uvel_bl.^2+vvel_bl.^2);
%        vel_bt_tmp = sqrt(uvel_bt.^2+vvel_bt.^2);
%        %vel_bl_tmp = abs(uvel_bl+vvel_bl)/2;
%        %vel_bt_tmp = abs(uvel_bt+vvel_bt)/2;
% 
%        velratio1 = abs((vel_bl-vel_bt)./vel_bl);
%        velratio2 = abs((vel_bl_tmp-vel_bt_tmp)./vel_bl_tmp);
% 
%        uvelratio2 = abs((uvel_bl-uvel_bt)./uvel_bl);
%        vvelratio2 = abs((vvel_bl-vvel_bt)./vvel_bl);
%         uvelratio = abs(uvel_shear./uvel_bl);
%         vvelratio = abs(vvel_shear./vvel_bl);
%        velratio2 = abs((vel_bl-vel_bt)./vel_bl);
% 
%        bgradx = abs(bgradx);
%        bgrady = abs(bgrady);


% 
%        %grad = sqrt(bgradx.^2+bgrady.^2);
        bgrad = (gradbx_abs+gradby_abs)/2;

        figure(999); hold on; plot(bgrad,epsU,'o',"Color",color(gg,:));
% 
%        if gg==1
%            %grad = abs(bgrady);
%            % velratio = abs((uvel_bl-uvel_bt)./uvel_bl);
%        else
%            %grad = abs(bgradx);
%            %velratio =  abs((vvel_bl-vvel_bt)./vvel_bl);
%        end
%        
%         sqCd = sqrt(2.5e-3); %non-dim
%         GammT = 0.02;
%         epsT_movmean = movmean(epsT,5*12,'omitnan');
%         epsT2_movmean = @(x)(movmean(velratio_gradb.*x./(sqCd*GammT+velratio_gradb.*x),5*12,'omitnan'));

        %% Find optimal E0
%         options = optimoptions(@fmincon,...
%         "PlotFcn",{@optimplotx,@optimplotfval,@optimplotfirstorderopt});
%         J = @(E)(rmse(epsT2_movmean(E),epsT_movmean,"all","omitnan"));
%         E0 = 3.6e-2;
%         [x,Jval,exitflag,output] = fmincon(J,E0,-1,0,[],[],[],[],[],options);
%         E0 = x(1);
       
%        
%        I = find(~isnan(epsT));
% 
%        urhs = uvelratio.*E0.*grad./(sqCd*GammT+velratio2.*E0.*grad);
%        vrhs = vvelratio.*E0.*grad./(sqCd*GammT+velratio2.*E0.*grad);
%        %uvrhs = velratio.*E0.*bgradx./(sqCd*GammT+velratio.*E0.*bgradx);
%        %rhs2 = velratio2.*E0.*grad./(sqCd*GammT+velratio2.*E0.*grad);
% 
% 
        epsT_data=load("epsT_estimates.mat");

        if jj==1

            %%% FIGURE 1

            title(ax_fig2(gg),basinslegend(gg));

            yyaxis(ax_fig2(gg),"left");
            g(1)=plot(ax_fig2(gg),MITTime,movmean(epsT,5*12),'-k','linewidth',2);
            lim = ylim(ax_fig2(gg));
            g(2)=plot(ax_fig2(gg),MITTime,movmean(TminTf_bl/TminTf_DI(24),5*12),'-','color',[0.47 0.67 0.19],'linewidth',1);
            g(3)=plot(ax_fig2(gg),MITTime,movmean(TminTf_bl(24)./TminTf_DI,5*12),'-','color',[0.93 0.69 0.13],'linewidth',1);
            ylim(ax_fig2(gg),[0.55 0.85]);

            yyaxis(ax_fig2(gg),"right");
            g(4)=plot(ax_fig2(gg),MITTime,movmean(bgrad,5*12),'--','color',[0.47 0.67 0.19],'linewidth',1);
            ylim(ax_fig2(gg),[0.02 0.06]);

            %A=epsT;
            %B=velratio_gradb.*E0./(sqCd*GammT+velratio_gradb.*E0);
            %A(:)\B(:)
            %plot(ax_fig2(gg),epsT,velratio_gradb.*E0./(sqCd*GammT+velratio_gradb.*E0),'.r');
            %plot(ax_fig2(gg),[0.4 0.9],[0.4 0.9]);
%             MITTime_epsT = epsT_data.epsT.(basin).MITTime;
%             epsT_full = epsT_data.epsT.(basin).epsT;
%             epsT_movmean = movmean(epsT_full,5*12,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_vargrad_vareps;
%             epsTprime_vargrad_vareps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_vargrad_ceps;
%             epsTprime_vargrad_ceps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_cgrad_vareps;
%             epsTprime_cgrad_vareps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_cgrad_ceps;
%             epsTprime_cgrad_ceps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%     
%             %plot(ax_fig2(gg),MITTime_epsT,epsT_full,'color',[0.7 0.7 0.7],'linewidth',0.5); hold on;
%             %plot(ax_fig2(gg),MITTime_epsT,epsTprime_full,'-','color',[1 0.75 0.25],'linewidth',0.5);
%             yyaxis(ax_fig2(gg),"left");
%             h(1) = plot(ax_fig2(gg),MITTime_epsT,epsT_movmean,'color','k','linewidth',2); hold on;
%             h(2) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_vargrad_vareps_movmean,'-','color',[1 0.41 0.16],'linewidth',2);
%             h(4) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_cgrad_vareps_movmean,':','color',[1 0.45 0.74],'linewidth',2);
%             
% 
%             yyaxis(ax_fig2(gg),"right");
%             h(3) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_vargrad_ceps_movmean,':','color',[0 0.45 0.74],'linewidth',2);  
%             h(5) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_cgrad_ceps_movmean,':','color',[0 1 0.74],'linewidth',2);
%             
             grid(ax_fig2(gg),'on'); box(ax_fig2(gg),"on");
% 
             xlim(ax_fig2(gg),[datenum('01012010','ddmmyyyy') datenum('01012220','ddmmyyyy')]); 
             xticks(ax_fig2(gg),[datenum('01012015','ddmmyyyy') datenum('01012040','ddmmyyyy') datenum('01012065','ddmmyyyy')...
             datenum('01012090','ddmmyyyy') datenum('01012115','ddmmyyyy') datenum(['01012140'],'ddmmyyyy') ...
             datenum('01012165','ddmmyyyy') datenum('01012190','ddmmyyyy') datenum('01012215','ddmmyyyy')]);
%             
            if gg==1  
                
                xticklabels(ax_fig2(gg),{});
                yyaxis(ax_fig2(gg),"right"); yticklabels(ax_fig2(gg),{});
                %ylim(ax_fig2(gg),[0.65 0.825]);
                leg = legend(g(:),{'$\epsilon_T$','$\frac{T_{\star}}{T_{\star\rm{DI}}(t=0)}$',...
                    '$\frac{T_{\star}(t=0)}{T_{\star\rm{DI}}}$','$\nabla b$'},...
                    "Interpreter","latex",'Orientation', 'Horizontal','fontsize',12);
                leg.Layout.Tile = 'north';

            elseif gg==2

                xticklabels(ax_fig2(gg),{});
                yyaxis(ax_fig2(gg),"left"); yticklabels(ax_fig2(gg),{});
                yyaxis(ax_fig2(gg),"right"); l=ylabel(ax_fig2(gg),'$\nabla b$','Interpreter','latex');
                l.Position = [825440 0.017 -1.0000];
                %ylim(ax_fig2(gg),[0.68 0.78]);

            elseif gg==3 

                xticklabels(ax_fig2(gg),{'0','','50','','100','','150','','200'});
                yyaxis(ax_fig2(gg),"right"); yticklabels(ax_fig2(gg),{});
                %ylim(ax_fig2(gg),[0.58 0.81]);
                
            elseif gg==4

                xticklabels(ax_fig2(gg),{'0','','50','','100','','150','','200'});
                yyaxis(ax_fig2(gg),"left"); yticklabels(ax_fig2(gg),{});
                %ylim(ax_fig2(gg),[0.53 0.85]);
    
    
            end
% 
%             xlimits = xlim(ax_fig2(gg)); ylimits = ylim(ax_fig2(gg));
%             E0 = epsT_data.epsT.(basin).E0_optimal;
%             text(ax_fig2(gg),xlimits(2)-5*360,ylimits(1)+(ylimits(2)-ylimits(1))/25,"E0 = "+extractBefore(string(E0),6),'horizontalalignment','right','verticalalignment','bottom','margin',2,'BackgroundColor','w');
             %yyaxis(tlo_fig2,"left"); 
             ylabel(tlo_fig2,'$\epsilon_T$','interpreter','latex');
             %yyaxis(tlo_fig2,"right"); ylabel(tlo_fig2,'$\nabla b$','interpreter','latex');
             xlabel(tlo_fig2,'Time [a]');


             %%% FIGURE 3

             title(ax_fig3(gg),basinslegend(gg));

            g(1)=plot(ax_fig3(gg),MITTime,movmean(epsU,5*12),'-k','linewidth',2);
            lim = ylim(ax_fig3(gg));
            g(2)=plot(ax_fig3(gg),MITTime,movmean(Ustar_bl/TminTf_DI(24),5*12),'-','color',[0.47 0.67 0.19],'linewidth',1);
            g(3)=plot(ax_fig3(gg),MITTime,movmean(Ustar_bl(24)./TminTf_DI,5*12),'-','color',[0.93 0.69 0.13],'linewidth',1);
            %ylim(ax_fig3(gg),[0.55 0.85]);

            %A=epsT;
            %B=velratio_gradb.*E0./(sqCd*GammT+velratio_gradb.*E0);
            %A(:)\B(:)
            %plot(ax_fig2(gg),epsT,velratio_gradb.*E0./(sqCd*GammT+velratio_gradb.*E0),'.r');
            %plot(ax_fig2(gg),[0.4 0.9],[0.4 0.9]);
%             MITTime_epsT = epsT_data.epsT.(basin).MITTime;
%             epsT_full = epsT_data.epsT.(basin).epsT;
%             epsT_movmean = movmean(epsT_full,5*12,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_vargrad_vareps;
%             epsTprime_vargrad_vareps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_vargrad_ceps;
%             epsTprime_vargrad_ceps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_cgrad_vareps;
%             epsTprime_cgrad_vareps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%             epsTprime_full = epsT_data.epsT.(basin).epsTprime_cgrad_ceps;
%             epsTprime_cgrad_ceps_movmean = movmean(epsTprime_full,12*5,"omitnan");
%     
%             %plot(ax_fig2(gg),MITTime_epsT,epsT_full,'color',[0.7 0.7 0.7],'linewidth',0.5); hold on;
%             %plot(ax_fig2(gg),MITTime_epsT,epsTprime_full,'-','color',[1 0.75 0.25],'linewidth',0.5);
%             yyaxis(ax_fig2(gg),"left");
%             h(1) = plot(ax_fig2(gg),MITTime_epsT,epsT_movmean,'color','k','linewidth',2); hold on;
%             h(2) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_vargrad_vareps_movmean,'-','color',[1 0.41 0.16],'linewidth',2);
%             h(4) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_cgrad_vareps_movmean,':','color',[1 0.45 0.74],'linewidth',2);
%             
% 
%             yyaxis(ax_fig2(gg),"right");
%             h(3) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_vargrad_ceps_movmean,':','color',[0 0.45 0.74],'linewidth',2);  
%             h(5) = plot(ax_fig2(gg),MITTime_epsT,epsTprime_cgrad_ceps_movmean,':','color',[0 1 0.74],'linewidth',2);
%             
             grid(ax_fig3(gg),'on'); box(ax_fig3(gg),"on");
% 
             xlim(ax_fig3(gg),[datenum('01012010','ddmmyyyy') datenum('01012220','ddmmyyyy')]); 
             xticks(ax_fig3(gg),[datenum('01012015','ddmmyyyy') datenum('01012040','ddmmyyyy') datenum('01012065','ddmmyyyy')...
             datenum('01012090','ddmmyyyy') datenum('01012115','ddmmyyyy') datenum(['01012140'],'ddmmyyyy') ...
             datenum('01012165','ddmmyyyy') datenum('01012190','ddmmyyyy') datenum('01012215','ddmmyyyy')]);
%             
            if gg==1  
                
                xticklabels(ax_fig3(gg),{});
                yyaxis(ax_fig3(gg),"right"); yticklabels(ax_fig3(gg),{});
                %ylim(ax_fig2(gg),[0.65 0.825]);
                leg = legend(g(:),{'$\epsilon_U$','$\frac{U_{\star}}{T_{\star\rm{DI}}(t=0)}$',...
                    '$\frac{U_{\star}(t=0)}{T_{\star\rm{DI}}}$'},...
                    "Interpreter","latex",'Orientation', 'Horizontal','fontsize',12);
                leg.Layout.Tile = 'north';

            elseif gg==2

                xticklabels(ax_fig3(gg),{});
                %yyaxis(ax_fig3(gg),"left"); yticklabels(ax_fig3(gg),{});
                
                %ylim(ax_fig2(gg),[0.68 0.78]);

            elseif gg==3 

                xticklabels(ax_fig3(gg),{'0','','50','','100','','150','','200'});
                %ylim(ax_fig2(gg),[0.58 0.81]);
                
            elseif gg==4

                xticklabels(ax_fig3(gg),{'0','','50','','100','','150','','200'});
                %yyaxis(ax_fig3(gg),"left"); yticklabels(ax_fig3(gg),{});
                %ylim(ax_fig2(gg),[0.53 0.85]);
    
    
            end
% 
%             xlimits = xlim(ax_fig2(gg)); ylimits = ylim(ax_fig2(gg));
%             E0 = epsT_data.epsT.(basin).E0_optimal;
%             text(ax_fig2(gg),xlimits(2)-5*360,ylimits(1)+(ylimits(2)-ylimits(1))/25,"E0 = "+extractBefore(string(E0),6),'horizontalalignment','right','verticalalignment','bottom','margin',2,'BackgroundColor','w');
             %yyaxis(tlo_fig2,"left"); 
             ylabel(tlo_fig3,'$\epsilon_U \; [\rm{m\,s}^{-1}\rm{K}^{-1}]$','interpreter','latex');
             %yyaxis(tlo_fig2,"right"); ylabel(tlo_fig2,'$\nabla b$','interpreter','latex');
             xlabel(tlo_fig3,'Time [a]');
             
        end

        %ylim(ylimits_left);



%         figure(22222+jj);
% % 
%         subplot(2,2,gg); hold on; title(basin);
%         
%         yyaxis left;
% 
%         plot(MITTime,epsT_movmean,'color','k','linewidth',2);
%      
%         ylimits_left = ylim;
% 
%         yyaxis right
% 
%         epsT2_movmean = movmean(epsT_new,5*12,'omitnan');
%         plot(MITTime,epsT2_movmean,'-','color','r');
% 
%         %ylim(ylimits_left);
% 
%         grid on; box on;


%             if gg==1
%     
%                 
%     
%                 %epsT_cgrad = velratio2.*E0.*bgrad./(sqCd*GammT+velratio2.*E0.*bgrad);
%     
%     
%                 %epsT_new = vvelratio.*E0.*bgrady./(sqCd*GammT+vvelratio.*E0.*bgrady);
%                 %epsT_cgrad = vvelratio.*E0.*bgrady(1)./(sqCd*GammT+vvelratio.*E0.*bgrady(1));
%                 %epsT_cvel = vvelratio(1).*E0.*bgrady./(sqCd*GammT+vvelratio(1).*E0.*bgrady);
%                 
%                 
%                 %plot(MITTime,movmean(epsT_cgrad,5*12,'omitnan'),'--','color','b');
%                 %plot(MITTime,movmean(epsT_cvel,5*12,'omitnan'),'--','color','g');
%                 %plot(MITTime,movmean(epsT_new2,15*12,'omitnan'),'-','color','b');
%     
%             else
%     
%                 %epsT_new = velratio_gradb.*E0(ee)./(sqCd*GammT+velratio_gradb.*E0(ee));
%                 %epsT_new2 = velratio2.*E0.*bgrad./(sqCd*GammT+velratio2.*E0.*bgrad);
%     
%     %            epsT_new = uvelratio.*E0.*bgradx./(sqCd*GammT+uvelratio.*E0.*bgradx);
%     %            epsT_cgrad = uvelratio.*E0.*bgradx(1)./(sqCd*GammT+uvelratio.*E0.*bgradx(1));
%     %            epsT_cvel = uvelratio(1).*E0.*bgradx./(sqCd*GammT+uvelratio(1).*E0.*bgradx);
%                 
%                 %plot(MITTime,movmean(epsT_new,5*12,'omitnan'),'-','color',CME0(ee,:));
%     %            plot(MITTime,movmean(epsT_cgrad,5*12,'omitnan'),'-','color','b');
%     %            plot(MITTime,movmean(epsT_cvel,5*12,'omitnan'),'-','color','g');
%                 %plot(MITTime,movmean(epsT_new2,5*12,'omitnan'),'-','color','b');
% 
%             end

      
        
        
% 
% 
%         figure(1111111+jj);
% 
%         subplot(2,2,gg); hold on; title(basin);
% 
%         yyaxis left;
% 
%         %plot(epsT,bgrady,'ok');
%         %plot(epsT,bgradx,'.g');
%         
% %        yyaxis left;
% 
% % 
%  %       plot(MITTime,epsT./epsT(24),'color','k','linewidth',2);
% % 
% %        ylim([0.55 0.9]);
% % 
% %         
% % 
% % %        movmean_epsT = movmean(epsT,5*12,'omitnan');
% % %        
%    %      yyaxis right;
% % %        if gg==1
% % %            rhs = urhs;
% % %        else
% % %            rhs = vrhs;
% % %        end
% % % 
%         plot(MITTime,bgradx,'o','color','m');
%         plot(MITTime,bgrady,'o','color','b');
%         plot(MITTime,bgrad,'o','color','c');
% 
%         yyaxis right;
        %plot(MITTime,abs(uvel_bl-uvel_bt)./abs(uvel_bl),'-','color','g');
        %plot(MITTime,abs(vvel_bl-vvel_bt)./abs(vvel_bl),'-','color','c');
        %plot(MITTime,abs(uvel_bt),'.','color','g');
        %plot(MITTime,abs(vvel_bt),'.','color','c');
 %%%%       plot(MITTime,uvelratio./uvelratio(24),'--','color','g');
 %%%%       plot(MITTime,vvelratio./vvelratio(24),'--','color','c');


%        plot(MITTime,uvelratio./uvelratio(1),':','color','m');
%        plot(MITTime,vvelratio./vvelratio(1),':','color','b');

        %plot(MITTime,grad,'-','color','g','linewidth',2);

  %      ylim([0 0.055]);
%        movmean_rhs = movmean(rhs,5*12,'omitnan');
% 
%         epsT_tmp = epsT(:);
%         rhs_tmp = rhs(:);
%        [R,P] = corr(epsT_tmp(I),rhs_tmp(I));
%        plot(R*epsT_tmp(I),rhs_tmp(I),'o','color',color(gg,:));
% 
%        [R,P] = corr(movmean_epsT(:),movmean_rhs(:));
%        plot(R*movmean_epsT,movmean_rhs,'ok');
% 
%        xlim([0 0.9]); ylim([0 0.9]);

%        yyaxis left;
%        plot(MITTime,movmean_epsT,'color','k','linewidth',2);
% 
%        yyaxis right;
%        plot(MITTime,movmean_rhs,'color','k','linewidth',2);
% 
%        ylabel('uv rhs');



%        subplot(2,2,2); hold on;
% 
%        yyaxis left;
%        plot(MITTime,epsT,'color',color(gg,:));
% 
%        yyaxis right;
%        plot(MITTime,uvel_shear,'-','color','k');
%        plot(MITTime,vvel_shear,'-','color','m');
%        plot(MITTime,vel_shear,'-','color','k');
%        ylabel('vel_shear')
%        
%        subplot(2,2,3); hold on;
%        yyaxis left;
%        plot(MITTime,epsT,'color',color(gg,:));
% 
%        yyaxis right;
%        plot(MITTime,bgradx,'-','color','m');
%        plot(MITTime,bgrady,'-','color','k');
%        ylabel('grad');
% 
%        subplot(2,2,4); hold on;
% 
%        yyaxis left;
%        plot(MITTime,epsT,'color',color(gg,:));
% 
%        yyaxis right;
%        plot(MITTime,uvelratio,'-','color','k');
%        plot(MITTime,vvelratio,'--','color','m');
%        plot(MITTime,velratio,'-','color','b');
%        ylabel('velratio');
% 
%        [R,P] = corrcoef(epsT(I),urhs(I));
%        basin
%        P(1,2)
%        R(1,2)
%        [R,P] = corrcoef(epsT(I),rhs2(I));
%        P(1,2)
%        R(1,2)
       %plot(MITTime,TminTf,'--','color',color(gg,:));
       %plot(MITTime,grad,'-k');

       %yyaxis right; 
       %plot(MITTime,grad,'-k');
       

       %plot(MITTime,velratio.*E0.*grad./(sqCd*GammT+velratio.*E0.*grad),'-k');
%       laggedcorr = xcorr(TminTf,TminTf_IF,10*12,'normalized');
%       figure(111); hold on; plot([-10*12:10*12],laggedcorr,'-o','color',color(gg,:),'linewidth',1);
%       grid on; box on;
       
%        laggedcorr = xcorr(movmean(TminTf_IF,5*12),movmean(mu,5*12),10*12);
%        figure(333); hold on;
%        plot(MITTime_IF,movmean(TminTf_IF,5*12),'-','color',color(gg,:));
%        plot(MITTime,movmean(TminTf,5*12),'--','color',color(gg,:));
%        figure(444); hold on;
%        plot([-10*12:10*12],laggedcorr,'-o','color',color(gg,:),'linewidth',1);
%        grid on; box on;

        %plot(ax_fig(jj),MITTime_IF,TminTf_IF,'-','color',color(gg,:),'linewidth',1);%
        %plot(ax_fig(jj),MITTime,TminTf,'--','color',color(gg,:),'linewidth',0.5);%
        %plot(ax_fig((jj-1)*3+1),MITTime,mu.^2,'-','color',color(gg,:),'linewidth',1);

 %%min(color(gg,:)+0.6,[1 1 1])

        if gg==1
             patch(ax_fig(jj),"XData",[MITTime_IF_001(1)-1e4 MITTime_IF_001(1)-1e4 MITTime_IF(end)+1e4 MITTime_IF(end)+1e4 MITTime_IF_001(1)-1e4],...
                 "YData",[1 5 5 1 1],"FaceColor",[0.5 0.5 0.5],"LineStyle","none","FaceAlpha",0.2);
             patch(ax_fig(jj+2),"XData",[MITTime_IF_001(1)-1e4 MITTime_IF_001(1)-1e4 MITTime_IF(end)+1e4 MITTime_IF(end)+1e4 MITTime_IF_001(1)-1e4],...
                 "YData",[1 5 5 1 1],"FaceColor",[0.5 0.5 0.5],"LineStyle","none","FaceAlpha",0.2);
             patch(ax_fig(jj+4),"XData",[MITTime_IF_001(1)-1e4 MITTime_IF_001(1)-1e4 MITTime_IF(end)+1e4 MITTime_IF(end)+1e4 MITTime_IF_001(1)-1e4],...
                 "YData",[1 5 5 1 1],"FaceColor",[0.5 0.5 0.5],"LineStyle","none","FaceAlpha",0.2);
             patch(ax_fig(jj+6),"XData",[MITTime_IF_001(1)-1e4 MITTime_IF_001(1)-1e4 MITTime_IF(end)+1e4 MITTime_IF(end)+1e4 MITTime_IF_001(1)-1e4],...
                 "YData",[1 5 5 1 1],"FaceColor",[0.5 0.5 0.5],"LineStyle","none","FaceAlpha",0.2);
             patch(ax_fig(jj+8),"XData",[MITTime_IF_001(1)-1e4 MITTime_IF_001(1)-1e4 MITTime_IF(end)+1e4 MITTime_IF(end)+1e4 MITTime_IF_001(1)-1e4],...
                 "YData",[1 5 5 1 1],"FaceColor",[0.5 0.5 0.5],"LineStyle","none","FaceAlpha",0.2);
        end

         movmean_TminTf = movmean(TminTf_IF.^2,5*12,'omitnan');
         scale1 = movmean_TminTf(1);
         plot(ax_fig(jj+2),MITTime_IF_001,TminTf_IF_001.^2/scale1,'-','color',colorsoft(gg,:),'linewidth',1);
         plot(ax_fig(jj+2),MITTime_IF,TminTf_IF.^2/scale1,'-','color',colorsoft(gg,:),'linewidth',1);
         g1((jj-1)*4+gg)=plot(ax_fig(jj+2),MITTime_IF,movmean_TminTf/scale1,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
%              
%         movmean_TminTf = movmean(TminTf.^2,5*12,'omitnan');
%         scale11 = movmean_TminTf(1);
%         %plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime_IF_001,TminTf_IF_001.^2/scale1,'-','color',min(color(gg,:)+0.6,[1 1 1]),'linewidth',1);
%         plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime,TminTf.^2/scale11,'-','color',min(color(gg,:)+0.5,[1 1 1]),'linewidth',1);
%         g1((jj-1)*4+gg)=plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime,movmean_TminTf/scale11,'linestyle',runIDlinestyle{jj},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
  
        movmean_mu = movmean(mu.^2,5*12,'omitnan');
        scale2 = movmean_mu(1);
        plot(ax_fig(jj+4),MITTime_001,mu_001.^2/scale2,'-','color',colorsoft(gg,:),'linewidth',1);
        plot(ax_fig(jj+4),MITTime,mu.^2/scale2,'-','color',colorsoft(gg,:),'linewidth',1);
        g2((jj-1)*4+gg)=plot(ax_fig(jj+4),MITTime,movmean_mu/scale2,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj));      

        movmean_epsT = movmean(epsT,5*12,'omitnan');
        scale3 = movmean_epsT(1);
        plot(ax_fig(jj+6),MITTime_001,epsT_001/scale3,'-','color',colorsoft(gg,:),'linewidth',1);
        plot(ax_fig(jj+6),MITTime,epsT/scale3,'-','color',colorsoft(gg,:),'linewidth',1);
        g3((jj-1)*4+gg)=plot(ax_fig(jj+6),MITTime,movmean_epsT/scale3,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
       
        movmean_epsU = movmean(epsU,5*12,'omitnan');
        scale4 = movmean_epsU(1);
        plot(ax_fig(jj+8),MITTime_001,epsU_001/scale4,'-','color',colorsoft(gg,:),'linewidth',1);
        plot(ax_fig(jj+8),MITTime,epsU/scale4,'-','color',colorsoft(gg,:),'linewidth',1);
        g4((jj-1)*4+gg)=plot(ax_fig(jj+8),MITTime,movmean_epsU/scale4,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj));       
       
        movmean_melt = movmean(melt,5*12,'omitnan');
        scale5 = scale1*scale2*scale3*scale4*m0*(365.25*24*60*60);
        plot(ax_fig(jj),MITTime_001,melt_001/scale5,'-','color',colorsoft(gg,:),'linewidth',1);
        plot(ax_fig(jj),MITTime,melt/scale5,'-','color',colorsoft(gg,:),'linewidth',1);
        g5((jj-1)*4+gg)=plot(ax_fig(jj),MITTime,movmean_melt/scale5,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 

%         figure(222+jj); subplot(2,2,gg); hold on; plot((TminTf_IF.^2-movmean_TminTf)/scale1,(mu.^2-movmean_mu)/scale2,'o','color',color(gg,:));
%         figure(444+jj); subplot(2,2,gg); hold on; plot((movmean_TminTf)/scale1,(movmean_mu)/scale2,'o','color',color(gg,:));
% 
        [R,P] = corrcoef(TminTf_IF.^2-movmean(TminTf_IF.^2,5*12,'omitnan'),...
            mu.^2-movmean(mu.^2,5*12,'omitnan'),'rows','complete');
        disp("corrcoef(TminTf_IF^2,mu^2)+ "+runID(jj)+" "+basin+" high freq. "+R(1,2));
        [R,P] = corrcoef(movmean(TminTf_IF.^2,5*12,'omitnan'),movmean(mu.^2,5*12,'omitnan'),'rows','complete');
        disp("corrcoef(TminTf_IF^2,mu^2)+ "+runID(jj)+" "+basin+" low freq. "+R(1,2));

        [R,P] = corrcoef(TminTf_IF-movmean(TminTf_IF,5*12,'omitnan'),...
            TminTf_DI-movmean(TminTf_DI,5*12,'omitnan'),'rows','complete');
        disp("corrcoef(TminTf_IF,TminTf_DI)+ "+runID(jj)+" "+basin+" high freq. "+R(1,2));
        [R,P] = corrcoef(movmean(TminTf_IF,5*12,'omitnan'),movmean(TminTf_DI,5*12,'omitnan'),'rows','complete');
        disp("corrcoef(TminTf_IF,TminTf_DI)+ "+runID(jj)+" "+basin+" low freq. "+R(1,2));

        if jj==1
            figure(1111);
            subplot(2,2,gg); hold on;
            %plot(MITTime,TminTf_IF-movmean(TminTf_IF,5*12,'omitnan'),'color','k');
            %plot(MITTime,TminTf_DI-movmean(TminTf_DI,5*12,'omitnan'),'color','r');
            plot(MITTime,movmean(TminTf_IF,5*12,'omitnan'),'linewidth',2,'color','k');
            plot(MITTime,movmean(TminTf_DI,5*12,'omitnan'),'linewidth',2,'color','r');
        end
%               
         Time = MITTime; n = round(log(numel(MITTime))/log(2));
         dt = (MITTime(end)-MITTime(1))/2^n;
         Time = linspace(MITTime(1),MITTime(end),2^n);
% 
         I = find(~isnan(TminTf_IF)); 
         TminTf_IF_psd = interp1(MITTime(I),TminTf_IF(I),Time);
         I = find(~isnan(mu)); 
         mu_psd = interp1(MITTime(I),mu(I),Time);
         I = find(~isnan(epsU));
         epsU_psd = interp1(MITTime(I),epsU(I),Time);
         I = find(~isnan(melt));
         melt_psd = interp1(MITTime(I),melt(I),Time);
% 
         fs = 1/(dt/365.25);
%         Y=fft(TminTf_IF_psd.^2/scale1);
%         P2 = abs(Y/numel(Time));
%         P1 = P2(1:numel(Time)/2+1);
%         P1(2:end-1) = 2*P1(2:end-1);
%         f = fs*(0:numel(Time)/2)/numel(Time);
%         figure(555); hold on; plot(f,P1,'color',color(gg,:));
%         
        [pxy,f] = cpsd(TminTf_IF_psd.^2/scale1,melt_psd/scale5,[],[],[],fs);
        [pxx,f] = pwelch(TminTf_IF_psd.^2/scale1,[],[],[],fs);
        figure(333); subplot(2,2,gg); hold on; 
        psd_color = [0.85 0.33 0.1;0 0.45 0.74];
        %yyaxis left;
        plot(1./f,abs(pxy)./abs(pxx),'color',psd_color(jj,:)); 
        ylim([0 2]);
%         yyaxis right;
%         plot(1./f,abs(pxx),'--','color',psd_color(jj,:));
%         ylim([0 6*1e-3]);
        title(basin);
        xlim([60/365.25 25]);
        set(gca,'XScale','log');
        grid on; box on;
        xlabel('Period [yr]'); ylabel('Transfer amplitude');



%         if gg==1
%             figure(111);
%             subplot(1,2,jj); hold on; 
%             plot(movmean_TminTf/scale1,movmean_melt/scale5,'.','color',color(gg,:));
%             m = (movmean_TminTf(:)/scale1)\(movmean_melt(:)/scale5);
%             plot([0.5 1.8],m*[0.5 1.8],'color',color(gg,:));
%             plot(movmean_mu/scale2,movmean_melt/scale5,'.','color',color(gg,:));
%             m = (movmean_mu(:)/scale2)\(movmean_melt(:)/scale5);
%             plot([0.5 1.8],m*[0.5 1.8],'color',color(gg,:));
%             plot(movmean_epsT/scale3,movmean_melt/scale5,'.','color',color(gg,:));
%             m = (movmean_epsT(:)/scale3)\(movmean_melt(:)/scale5);
%             plot([0.5 1.8],m*[0.5 1.8],'color',color(gg,:));
%             plot(movmean_epsU/scale4,movmean_melt/scale5,'.','color',color(gg,:));
%             m = (movmean_epsU(:)/scale4)\(movmean_melt(:)/scale5);
%             plot([0.5 1.8],m*[0.5 1.8],'color',color(gg,:));
%             plot([0 2],[0 2],'-k');
%             ylim([0.5 1.8]);
%             xlim([0.5 1.8]);
%             %plot(MITTime,m0*(365.25*24*60*60)*scale1*scale2*epsT.*epsU/scale5);
%             %plot(MITTime,m0*(365.25*24*60*60)*scale2*TminTf_IF.^2.*epsT.*epsU/scale5);
%             %plot(MITTime,m0*(365.25*24*60*60)*scale3*TminTf_IF.^2.*mu.^2.*epsU/scale5);
%             %plot(MITTime,m0*(365.25*24*60*60)*scale4*TminTf_IF.^2.*mu.^2.*epsT/scale5);
%             %plot(MITTime,movmean_melt/scale5,'linewidth',2);
%         end
%         
%         MAPE_Tstar2(gg,1) = mape(m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_mu2(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_epsT(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_epsU(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4/scale5,melt/scale5,'omitnan');

        RMSE_Tstar2(gg,1) = rmse(m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
        RMSE_mu2(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
        RMSE_epsT(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU/scale5,melt/scale5,'omitnan');
        RMSE_epsU(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4/scale5,melt/scale5,'omitnan');

         figure(200+jj*10+gg); hold on;
         t1 = melt;
         t2 = m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU;
         plot(MITTime,(t2-t1)./t2);
         t2 = m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU;
         plot(MITTime,(t2-t1)./t2);
         t2 = m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU;
         plot(MITTime,(t2-t1)./t2);
         t2 = m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4;
         plot(MITTime,(t2-t1)./t2);
         %plot(MITTime,melt);
         title([basin,' - ',runID{jj}]);


%         [R,P] = corrcoef(movmean_TminTf,movmean_melt);
%         R2_Tstar2(gg,1) = R(1,2)^2*100;
%         P_Tstar2(gg,1) = P(1,2);
         scale_Tstar(gg,1) = sqrt(scale1);
%         [R,P] = corrcoef(movmean_mu,movmean_melt);
%         R2_mu2(gg,1) = R(1,2)^2*100;
%         P_mu2(gg,1) = P(1,2);
         scale_mu(gg,1) = sqrt(scale2);
%         [R,P] = corrcoef(movmean_epsU,movmean_melt);
%         R2_epsU(gg,1) = R(1,2)^2*100;
%         P_epsU(gg,1) = P(1,2);
         scale_epsU(gg,1) = scale4;
%         [R,P] = corrcoef(movmean_epsT,movmean_melt);
%         R2_epsT(gg,1) = R(1,2)^2*100;
%         P_epsT(gg,1) = P(1,2);
         scale_epsT(gg,1) = scale3;
         scale_melt(gg,1) = scale5;


%         g6(gg)=plot(ax_fig2(floor(gg/3)+1),TminTf,melt,'o','color',color(gg,:),'markersize',2,'MarkerFaceColor',color(gg,:));
%         g7(gg)=plot(ax_fig2(floor(gg/3)+1),TminTf_IF,melt,'d','color',color(gg,:),'markersize',4,'MarkerFaceColor',min(color(gg,:)+0.5,[1 1 1]),'MarkerEdgeColor','none');
%             TminTf_av = [3 4 3 2;6 6 4 5];
%             figure(111); hold on;
%             subplot(1,4,gg); hold on;
%             plot(MITTime,melt,'-','color',min(color(gg,:)+0.6-(jj-1)/5,[1 1 1]),'linewidth',1);
%             plot(MITTime,m0*mu.^2.*epsT.*epsU.*TminTf_av(jj,gg)*(365.25*24*60*60),'-','color',color(gg,:),'linewidth',1.5);
%         g9 = plot(ax_fig2(1),-10,-10,'o','color',[0 0 0],'markersize',3,'MarkerFaceColor',[0 0 0]);
%         g10 = plot(ax_fig2(1),-10,-10,'o','color',[0.7 0.7 0.7],'markersize',6,'MarkerFaceColor',[0.7 0.7 0.7]);
%         g11 = plot(ax_fig2(1),-10,-10,'-','color',[0 0 0],'linewidth',1.5);
% 
%         g12(gg)=plot(ax_fig3(2),MITTime,vvel,'color',color(gg,:),'linewidth',1.5);
%         plot(ax_fig3(3),MITTime,uvel,'color',color(gg,:),'linewidth',1.5);
%         plot(ax_fig3(4),MITTime,bgradx,'color',color(gg,:),'linewidth',1.5);
%         plot(ax_fig3(5),MITTime,bgrady,'color',color(gg,:),'linewidth',1.5);
%         %plot(ax_fig3(6),MITTime,gradxrho,'color',color(gg,:),'linewidth',1.5);
%         %plot(ax_fig3(7),MITTime,gradyrho,'color',color(gg,:),'linewidth',1.5);
%         plot(ax_fig3(6),MITTime,overturning_IF/1e6,'color',color(gg,:)); % in Sv
% 
%         nt = numel(epsU);
%         facecolor = [linspace(1,color(gg,1),ceil(1.3*nt))',...
%             linspace(1,color(gg,2),ceil(1.3*nt))',linspace(1,color(gg,3),ceil(1.3*nt))'];
%         facecolor = facecolor([end-nt+1:end],:);
%         scatter(ax_fig4(1),bsf/1e6,epsU,13,facecolor,'filled');
%         scatter(ax_fig4(2),bsf/1e6,melt,13,facecolor,'filled');
%         scatter(ax_fig4(3),osf/1e6,melt,13,facecolor,'filled');

    end

    T = table(basins',scale_Tstar,RMSE_Tstar2,scale_mu,RMSE_mu2,scale_epsT,RMSE_epsT,scale_epsU,RMSE_epsU,scale_melt)


end

title(ax_fig(1),'$\mbox{\textit{hi\_melt}}$','fontsize',12,'interpreter','latex');
title(ax_fig(2),'$\mbox{\textit{av\_melt}}$','fontsize',12,'interpreter','latex');


for kk=1:5*numel(runID)

    patch(ax_fig(kk),[datenum('01011980','ddmmyyyy') datenum('01011980','ddmmyyyy') ...
        datenum('01012015','ddmmyyyy') datenum('01012015','ddmmyyyy')],[-10 100 100 -10],'w','facealpha',0,'linestyle','--','linewidth',2,'Edgecolor',[0.6 0.6 0.6]);
    xlim(ax_fig(kk),[datenum('01011990','ddmmyyyy') datenum('01012215','ddmmyyyy')]); 
    xticks(ax_fig(kk),[datenum('01012015','ddmmyyyy') datenum('01012040','ddmmyyyy') datenum('01012065','ddmmyyyy')...
        datenum('01012090','ddmmyyyy') datenum('01012115','ddmmyyyy') datenum(['01012140'],'ddmmyyyy') ...
        datenum('01012165','ddmmyyyy') datenum('01012190','ddmmyyyy') datenum('01012215','ddmmyyyy')]);

    h(1) = plot(ax_fig(kk),[0 1],[-1 -1],'-k','linewidth',2);
    h(2) = plot(ax_fig(kk),[0 1],[-1 -1],'--k','linewidth',2);

    text(ax_fig(kk),datenum('01011995','ddmmyyyy'),2.4,panellabel{kk},"HorizontalAlignment",'left','VerticalAlignment','top','fontsize',13);
    
    %datetick(ax_fig(kk),'x','keeplimits'); 
    if kk<9
        xticklabels(ax_fig(kk),{''});
    else
        xticklabels(ax_fig(kk),{'0','','50','','100','','150','','200'});
        %xlabel(ax_fig(kk),'Time [yrs]','FontSize',12);
    end
    if kk<3
        ylim(ax_fig(kk),[0 2.5]);
        yticks(ax_fig(kk),[0:0.25:2.5]);       
    end
    if ismember(kk,[3:4])
        ylim(ax_fig(kk),[0 2.5]);
        yticks(ax_fig(kk),[0:0.25:2.5]);
    end
    if ismember(kk,[5:6])
        ylim(ax_fig(kk),[0 2.5]);
        yticks(ax_fig(kk),[0:0.25:2.5])
    end
    if ismember(kk,[7:8])
        ylim(ax_fig(kk),[0 2.5]);
        yticks(ax_fig(kk),[0:0.25:2.5]);
    end
    if ismember(kk,[9:10])
        ylim(ax_fig(kk),[0 2.5]);
        yticks(ax_fig(kk),[0:0.25:2.5]);
    end
    if mod(kk,2)==1
        ylabel(ax_fig(kk),ylabels{(kk-1)/2+1},'interpreter','latex','fontsize',12);
        yticklabels(ax_fig(kk),{'0','','0.5','','1','','1.5','','2','','2.5'});
    else
        yticklabels(ax_fig(kk),{''});
    end
    grid(ax_fig(kk),'on'); box(ax_fig(kk),'on');    
end

%     for kk=1:2
%         xlabel(ax_fig2(kk),'$\overline{T_{\star}}\;\;\rm{[}{}^{\circ}\rm{C]}$','interpreter','latex');
%         xlim(ax_fig2(kk),[0 3]); ylim(ax_fig2(kk),[0 70]);
%         grid(ax_fig2(kk),'on'); box(ax_fig2(kk),'on');
%     end
%     yticklabels(ax_fig2(2),{''});
%     ylabel(ax_fig2(1),'$\overline{m}\;\;\rm{[m a}{}^{-1}\rm{]}$','interpreter','latex');

    
%     if ~contains(runID,'PTDC_004')
%         load(['/Volumes/mainJDRydt2/UaMITgcm_v2/cases/',runID{1},'/output/202001/Ua/',runID{1},'-RestartFile.mat']);
%     else
%         load('/Volumes/mainJDRydt2/UaMITgcm_v2/cases/PTDC_001/output/201501/Ua/PTDC_001-RestartFile.mat');
%     end
%     GF.node(GF.node<1)=0; %GF.node(GF.node==0 & F.b<-400)= -1;
% %     PlotNodalBasedQuantities_JDR(ax_fig3(1),MUA.connectivity,MUA.coordinates,...
% %         GF.node,CtrlVarInRestartFile,'EdgeColor','none'); 
% %     colormap(ax_fig3(1),[0.8 0.8 0.8 ; 1 1 1]);
% %     GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVarInRestartFile);
% %     xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
% %     [xGL,yGL]=LineUpEdges2([],xa,xb,ya,yb);
% %     
% %     if ~contains(runID,'PTDC_004')
% %         plot(ax_fig3(1),xGL,yGL,'-b','linewidth',2); 
% %         load(['/Volumes/mainJDRydt2/UaMITgcm_v2/cases/',runID{1},'/output/201501/Ua/',runID{1},'-RestartFile.mat']);
% %         GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVarInRestartFile);
% %         xa=GLgeo(:,3) ;  xb=GLgeo(:,4) ; ya=GLgeo(:,5) ;  yb=GLgeo(:,6) ;
% %         [xGL,yGL]=LineUpEdges2([],xa,xb,ya,yb);
% %     end
% 
%     %plot(ax_fig3(1),xGL,yGL,'color',[133 87 35]/255,'linewidth',2);
% %     plot(ax_fig3(1),MUA.Boundary.x,MUA.Boundary.y,'-k');
% %     dL = sqrt((section(1).xmid(2:end)-section(1).xmid(1:end-1)).^2+(section(1).ymid(2:end)-section(1).ymid(1:end-1)).^2);
% %     J=find(dL>2e3); section(1).xmid(J)=NaN; section(1).ymid(J)=NaN;
% %     plot(ax_fig3(1),section(1).xmid,section(1).ymid,'.-','color',[133 87 35]/255,'markersize',2,'markerfacecolor',[133 87 35]/255);
% 
%     if ~contains(runID,'PTDC_004')
%         [~,I] = min(abs(MITTime-datenum('01012500','ddmmyyyy')));
%         dL = sqrt((section(I).xmid(2:end)-section(I).xmid(1:end-1)).^2+(section(I).ymid(2:end)-section(I).ymid(1:end-1)).^2);
%         J=find(dL>2e3); section(I).xmid(J)=NaN; section(I).ymid(J)=NaN;
%         plot(ax_fig3(1),section(I).xmid,section(I).ymid,'.-b','markersize',2,'markerfacecolor','b');
%     end
% 
%     axis(ax_fig3(1),'equal');
%     xlim(ax_fig3(1),[ -1700e3 -1450e3]);
%     ylim(ax_fig3(1),[-700e3 -230e3]);
%     grid(ax_fig3(1),'off'); box(ax_fig3(1),'off'); yticklabels(ax_fig3(1),{''}); xticklabels(ax_fig3(1),{''});
% % pos = get(H2,'Position');
% set(H2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% fname = ['../Figures/',runID{1},'/Tsquare'];
% print(H2,fname,'-dpng','-r400');
% 
%     for kk=2:6
%         xlim(ax_fig3(kk),[datenum('01012015','ddmmyyyy') datenum('01012500','ddmmyyyy')]); 
%         datetick(ax_fig3(kk),'x','keeplimits'); 
%         grid(ax_fig3(kk),'on'); box(ax_fig3(kk),'on');   
%     end
%     for kk=2:5
%         xticklabels(ax_fig3(kk),{''});
%     end
%     ylabel(ax_fig3(2),{'$U_y\;\; \rm{[m/s]}$'},'interpreter','latex'); ylim(ax_fig3(2),[-0.07 0.12]);
%     ylabel(ax_fig3(3),{'$U_x\;\; \rm{[m/s]}$'},'interpreter','latex'); ylim(ax_fig3(3),[-0.07 0.12]);
%     ylabel(ax_fig3(4),{'$\partial_x b$'},'interpreter','latex'); ylim(ax_fig3(4),[-50 15]*1e-3);
%     ylabel(ax_fig3(5),{'$\partial_y b$'},'interpreter','latex'); ylim(ax_fig3(5),[-50 15]*1e-3);
%     %ylabel(ax_fig3(6),{'$\partial_x \rho$'},'interpreter','latex'); ylim(ax_fig3(6),[-15 5]*1e-6);
%     %ylabel(ax_fig3(7),{'$\partial_y \rho$'},'interpreter','latex'); ylim(ax_fig3(7),[-15 5]*1e-6);
%     ylabel(ax_fig3(6),{'Overturning [Sv]'});
%     title(ax_fig3(2),{'x-momentum'});
%     title(ax_fig3(3),{'y-momentum'});
%     title(ax_fig3(6),'Overturning');
% 
%     xlim(ax_fig4(1),[0 1.3]); ylim(ax_fig4(1),[0 3.5]*1e-3);
%     grid(ax_fig4(1),'on'); box(ax_fig4(1),'on');
%     for kk=2:3
%         xlim(ax_fig4(kk),[0 1.3]);
%         ylim(ax_fig4(kk),[0 85]);
%         grid(ax_fig4(kk),'on'); box(ax_fig4(kk),'on');
%     end    
%     xlabel(ax_fig4(1),'$\mbox{Barotropic streamfunction [Sv]}$','interpreter','latex');
%     ylabel(ax_fig4(1),'$\overline{\epsilon_U}(t)\;\;\rm{[m\,s}^{-1}{}^{\circ}\rm{C}^{-1}\rm{]}$','interpreter','latex');
%     xlabel(ax_fig4(2),'$\mbox{Barotropic streamfunction [Sv]}$','interpreter','latex');
%     ylabel(ax_fig4(2),'$\overline{m}\;\;\rm{[m/yr]}$','interpreter','latex');
%     xlabel(ax_fig4(3),'$\mbox{Overturning streamfunction [Sv]}$','interpreter','latex');
%     ylabel(ax_fig4(3),'$\overline{m}\;\;\rm{[m/yr]}$','interpreter','latex');

for jj=1:numel(runID)
    for gg=[1 3 4 2]
        uistack(g1((jj-1)*4+gg),'top');
        uistack(g2((jj-1)*4+gg),'top');
        uistack(g3((jj-1)*4+gg),'top');
        uistack(g4((jj-1)*4+gg),'top');
        uistack(g5((jj-1)*4+gg),'top');
%        uistack(g6(gg),'top');
%        uistack(g8(gg),'top');
    end
end

leg1=legend(g1(1:4),basinslegend{:},'Orientation', 'Horizontal','fontsize',10);
leg1.Layout.Tile = 'north';

%ah1=axes('position',get(gca,'position'),'visible','off');

%leg2=legend(ah1,h(1:2),{'$\mbox{\textit{hi\_melt}}$','$\mbox{\textit{av\_melt}}$'},...
%    'Orientation', 'Horizontal','fontsize',14,'interpreter','latex','FontWeight','bold');
%leg2.Position = [0.1 0.95 0.14 0.02];

xlabel(tlo_fig,'time [yr]','fontsize',10);
% 
% leg1=legend(ax_fig2(2),g8(:),basinslegend{:},'Orientation', 'Horizontal');
% leg1.Layout.Tile = 'north';
% 
% leg2=legend(ax_fig2(1),[g10 g9 g11],{'$\overline{T_{\star{\rm in,IF}}}$','$\overline{T_{\star{\rm in}}}$',...
%     '$m_0\,\eps_U_0\,\eps_T_0\left(\overline{T_{\star,{\rm in}}}\right)^2$'},'Orientation','Vertical','Location','northwest','Interpreter','latex');
% 
% leg=legend(g12(:),basinslegend{:},'Orientation', 'Horizontal');
% leg.Layout.Tile = 'north';
% 
% leg3=legend(ax_fig3(2),g12(:),basinslegend{:},'Orientation', 'Vertical','Location','northwest');

pos = get(H1,'Position');
set(H1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/Fig5_GeometricCoefficients'];
print(H1,fname,'-dpng','-r400');

pos = get(H2,'Position');
set(H2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/FigS3_epsT_gradb_',runID{1}];
print(H2,fname,'-dpng','-r400');

pos = get(H3,'Position');
set(H3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/FigS4_epsU_',runID{1}];
print(H3,fname,'-dpng','-r400');

% pos = get(H2,'Position');
% set(H2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% fname = ['../Figures/',runID{1},'/Tsquare'];
% print(H2,fname,'-dpng','-r400');
% 
% pos = get(H3,'Position');
% set(H3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% fname = ['../Figures/',runID{1},'/GeostrophicBalance'];
% print(H3,fname,'-dpng','-r400');
% 
% pos = get(H4,'Position');
% set(H4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
% fname = ['../Figures/',runID{1},'/Streamfunctions'];
% print(H4,fname,'-dpng','-r400');