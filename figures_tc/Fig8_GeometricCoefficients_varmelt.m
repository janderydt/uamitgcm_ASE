function Fig8_GeometricCoefficients_varmelt

addpath(getenv("froot_tools"));

endtime = datenum('01012215','ddmmyyyy');

runID = {'PTDC_000','PTDC_001','PTDC_003'};
runIDlinestyle = {'-','-','-'};
runIDlinewidth = [2 2 2];
basins = {'PIG','TW','CR','DT'};
basinslegend = {'Pine Island','Thwaites','Crosson','Dotson'};
%color = [[51 34 136];[17 119 51];[221 204 119];[136 34 85]]/255;
%color = [[51 34 136];[17 119 51];[133 87 35];[136 34 85]]/255;
%color=[[166,206,227];[31,120,180];[178,223,138];[51,160,44]]/255;
%color = [[217,95,2];[117,112,179];[231,41,138];[102,166,30]]/255;
%color = [[228,26,28];[55,126,184];[77,175,74];[152,78,163]]/255;
%colorsoft = [[251,180,174];[179,205,227];[204,235,197];[222,203,228]]/255;
%color = [[31,120,180];[27,158,119];[246,190,0];[152,78,163]]/255;
%colorsoft = [[166,206,227];[141,211,199];[251,219,101];[190,174,212]]/255;

%color = [166,206,227;31,120,180]/255;
%colorsoft = [128 128 128;51 51 51]/255;
color = [0.47 0.67 0.19;0.2 0.2 0.2];
colorsoft = [0.47 0.67 0.19;0.2 0.2 0.2;0.85 0.33 0.1];

rhoConst = 1024;
Cp = 3974;
Lf = 334e3; %J/kg
Cd = 0.0025;
GammT = 0.015;
m0 = 2.2e-4;

ylabels = {'$\widetilde{T}_{\star {\rm IF}}{}^2$',...
    '$\widetilde{\mu}^2$','$\widetilde{\epsilon}_T$','$\widetilde{\epsilon}_U$',...
    '$\widetilde{m}$'};

H1=fig('units','inches','width',70*12/72.27,'height',60*12/72.27,'fontsize',12,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);

tlo_fig = tiledlayout(2,2,'TileSpacing','tight');
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

% H2 = fig('units','inches','width',60*12/72.27,'height',35*12/72.27,'fontsize',16,'font','Helvetica','defaultAxesColorOrder',[0 0 0; 0 0 0]);
% 
% tlo_fig2 = tiledlayout(1,2,'TileSpacing','tight');
% for i = 1:2
%     ax_fig2(i) = nexttile(tlo_fig2); hold on;
% end
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
%     ax_fig4(i) = nexttile(tlo_fig4); hold on;
% end

for jj=1:numel(runID)

    load(['HeatVolumeTransport_IceFront_below400m_',runID{jj},'.mat']);
    MITTime_IF = MITTime;
    section_IF = section;

    load(['HeatVolumeTransport_moving400mdraft_',runID{jj},'.mat']);
       
    %SF = load('CavityStreamFunctions_PTDC_002_monthly.mat');

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

        TminTf_IF=[]; kk=1;


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

       Ustar_bl =[]; TminTf_bl=[]; TminTf=[]; melt=[]; 
       bgradx=[]; bgrady=[]; uvel=[]; vvel=[]; 
       gradxrho=[]; gradyrho=[]; uvel_bt=[]; vvel_bt=[];
       kk=1;

       for tt=start_MITTime:end_MITTime

            I = find(inpoly([section(tt).xmid(:) section(tt).ymid(:)],[xpoly(:) ypoly(:)]));
            
            % thermal driving, heat flux, volume flux, cross section
            %TminTf(kk) = mean(section(tt).monthly.TminTf_belowz_in_depthmean(I),'all','omitnan');
            TminTf(kk) = mean(section(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');

            Ustar_bl(kk) = integral2D.(basin).monthly.blUStar_mean(tt);
            TminTf_bl(kk) = integral2D.(basin).monthly.TminTf_mean(tt);

            melt(kk) = integral2D.(basin).monthly.melt_integral(tt)./integral2D.(basin).monthly.ISarea_integral(tt)*1e12/1e3; %m/yr

            %bgrady(kk)=integral2D.(basin).monthly.bgrady_mean(tt);
            %bgradx(kk)=integral2D.(basin).monthly.bgradx_mean(tt);
            %vvel(kk)=integral2D.(basin).monthly.blVVEL_mean(tt);
            %uvel(kk)=integral2D.(basin).monthly.blUVEL_mean(tt);
            
            %uvel_bt(kk) = integral2D.(basin).monthly.btUVEL_mean(tt);
            %vvel_bt(kk) = integral2D.(basin).monthly.btVVEL_mean(tt);
            %gradxrho(kk)=integral2D.(basin).monthly.gradxRho_mean(tt);
            %gradyrho(kk)=integral2D.(basin).monthly.gradyRho_mean(tt);

            kk= kk+1;
    
       end

       MITTime = MITTime(start_MITTime:end_MITTime);
       [MITTime,ia,ic] = unique(MITTime);
       TminTf = TminTf(ia); 
       melt = melt(ia);
       Ustar_bl = Ustar_bl(ia);
       TminTf_bl = TminTf_bl(ia);
       %bgradx = bgradx(ia);
       %bgrady = bgrady(ia);
       %uvel = uvel(ia);
       %vvel = vvel(ia);
       %gradxrho = gradxrho(ia);
       %gradyrho = gradyrho(ia);

       %figure(222+jj); hold on;
       %plot(MITTime,TminTf,'color',color(gg,:));
       %plot(MITTime_IF,TminTf_IF,'--','color',color(gg,:))

       mu = TminTf./interp1(MITTime_IF,TminTf_IF,MITTime);
       epsU = Ustar_bl./TminTf;
       epsT = TminTf_bl./TminTf;
%        bsf = interp1(SF.time,SF.bsf.(basin),MITTime);
%        osf = interp1(SF.time,SF.osf.(basin),MITTime);

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

         movmean_TminTf = movmean(TminTf_IF.^2,5*12,'omitnan');
         scale1 = movmean_TminTf(1);
         %plot(ax_fig(jj),MITTime_IF_001,TminTf_IF_001.^2,'-','color',colorsoft(gg,:),'linewidth',1);
         %plot(ax_fig(jj),MITTime_IF,TminTf_IF.^2,'-','color',colorsoft(gg,:),'linewidth',1);
         %g1((jj-1)*4+gg)=plot(ax_fig(gg),MITTime_IF,movmean_TminTf,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
%              
%         movmean_TminTf = movmean(TminTf.^2,5*12,'omitnan');
%         scale11 = movmean_TminTf(1);
%         %plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime_IF_001,TminTf_IF_001.^2/scale1,'-','color',min(color(gg,:)+0.6,[1 1 1]),'linewidth',1);
%         plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime,TminTf.^2/scale11,'-','color',min(color(gg,:)+0.5,[1 1 1]),'linewidth',1);
%         g1((jj-1)*4+gg)=plot(ax_fig(2*(jj-1)+floor(gg/5)+1),MITTime,movmean_TminTf/scale11,'linestyle',runIDlinestyle{jj},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
  
        movmean_mu = movmean(mu.^2,5*12,'omitnan');
        scale2 = movmean_mu(1);
        %plot(ax_fig(jj+2),MITTime_001,mu_001.^2,'-','color',colorsoft(gg,:),'linewidth',1);
        %plot(ax_fig(jj+2),MITTime,mu.^2,'-','color',colorsoft(gg,:),'linewidth',1);
        %g2((jj-1)*4+gg)=plot(ax_fig(gg),MITTime,movmean_mu,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj));      

        movmean_epsT = movmean(epsT,5*12,'omitnan');
        scale3 = movmean_epsT(1);
        %plot(ax_fig(jj+4),MITTime_001,epsT_001,'-','color',colorsoft(gg,:),'linewidth',1);
        %plot(ax_fig(jj+4),MITTime,epsT,'-','color',colorsoft(gg,:),'linewidth',1);
        %g3((jj-1)*4+gg)=plot(ax_fig(gg),MITTime,movmean_epsT,'linestyle',runIDlinestyle{gg},'color',color(gg,:),'linewidth',runIDlinewidth(jj)); 
       
        if contains(runID{jj},{'PTDC_000','PTDC_001'})
            %yyaxis(ax_fig(gg),"right")
            movmean_epsU = movmean(epsU,5*12,'omitnan');
            scale4 = movmean_epsU(1);
            %plot(ax_fig(jj+6),MITTime_001,epsU_001,'-','color',colorsoft(gg,:),'linewidth',1);
            %plot(ax_fig(jj+6),MITTime,epsU,'-','color',colorsoft(gg,:),'linewidth',1);
            %g1(jj)=plot(ax_fig(gg),MITTime,movmean_epsU,'linestyle',':','color',color(jj,:),'linewidth',runIDlinewidth(jj));   
        end
       
        %yyaxis(ax_fig(gg),"left")
        movmean_melt = movmean(melt,5*12,'omitnan');
        scale5 = scale1*scale2*scale3*scale4*m0*(365.25*24*60*60);
        %plot(ax_fig(jj+8),MITTime_001,melt_001,'-','color',colorsoft(gg,:),'linewidth',1);
        %plot(ax_fig(jj+8),MITTime,melt,'-','color',colorsoft(gg,:),'linewidth',1);
        g2(jj)=plot(ax_fig(gg),MITTime,movmean_melt,'linestyle','-','color',colorsoft(jj,:),'linewidth',runIDlinewidth(jj)); 

        meltdata(jj).(basin).movmean_melt = movmean_melt;
        meltdata(jj).(basin).time = MITTime;
        %figure(222+jj); subplot(2,2,gg); hold on; plot((TminTf_IF.^2-movmean_TminTf)/scale1,(mu.^2-movmean_mu)/scale2,'o','color',color(gg,:));
        %figure(444+jj); subplot(2,2,gg); hold on; plot((movmean_TminTf)/scale1,(movmean_mu)/scale2,'o','color',color(gg,:));

        [R,P] = corrcoef(TminTf_IF.^2-movmean_TminTf,mu.^2-movmean_mu,'rows','complete');
        disp(runID(jj)+" "+basin+" high freq. "+R(1,2));
        [R,P] = corrcoef(movmean_TminTf,movmean_mu,'rows','complete');
        disp(runID(jj)+" "+basin+" low freq. "+R(1,2));
              
        Time = MITTime; n = round(log(numel(MITTime))/log(2));
        dt = (MITTime(end)-MITTime(1))/2^n;
        Time = linspace(MITTime(1),MITTime(end),2^n);

        I = find(~isnan(TminTf_IF)); 
        TminTf_IF_psd = interp1(MITTime(I),TminTf_IF(I),Time);
        I = find(~isnan(mu)); 
        mu_psd = interp1(MITTime(I),mu(I),Time);
        I = find(~isnan(epsU));
        epsU_psd = interp1(MITTime(I),epsU(I),Time);
        I = find(~isnan(melt));
        melt_psd = interp1(MITTime(I),melt(I),Time);

%         fs = 1/(dt/365.25);
%         Y=fft(TminTf_IF_psd.^2/scale1);
%         P2 = abs(Y/numel(Time));
%         P1 = P2(1:numel(Time)/2+1);
%         P1(2:end-1) = 2*P1(2:end-1);
%         f = fs*(0:numel(Time)/2)/numel(Time);
% %        figure(555); hold on; plot(f,P1,'color',color(gg,:));
%         
%         [pxy,f] = cpsd(TminTf_IF_psd.^2/scale1,melt_psd/scale5,[],[],[],fs);
%         [pxx,f] = pwelch(TminTf_IF_psd.^2/scale1,[],[],[],fs);
%         figure(333+jj); subplot(2,2,gg); hold on; 
%         yyaxis left;
%         plot(1./f,abs(pxy)./abs(pxx),'color',color(gg,:)); 
%         %ylim([0 1.5]);
%         yyaxis right;
%         plot(1./f,abs(pxx),'--','color',color(gg,:));
%         %ylim([0 6*1e-3]);
%         title(basin);
%         xlim([60/365.25 25]);
%         set(gca,'XScale','log');
%         grid on; box on;
%         xlabel('Period (years)'); ylabel('Transfer amplitude');



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
        
%         MAPE_Tstar2(gg,1) = mape(m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_mu2(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_epsT(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU/scale5,melt/scale5,'omitnan');
%         MAPE_epsU(gg,1) = mape(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4/scale5,melt/scale5,'omitnan');
% 
%         RMSE_Tstar2(gg,1) = rmse(m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU,melt,'omitnan');
%         RMSE_mu2(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU,melt,'omitnan');
%         RMSE_epsT(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU,melt,'omitnan');
%         RMSE_epsU(gg,1) = rmse(m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4,melt,'omitnan');
% 
% %         figure(222+gg); hold on;
% %         plot(MITTime,m0*(365.25*24*60*60)*scale1*mu.^2.*epsT.*epsU/scale5,MITTime,melt/scale5);
% %         hold on; plot(MITTime,m0*(365.25*24*60*60)*TminTf_IF.^2.*scale2.*epsT.*epsU/scale5);
% %         hold on; plot(MITTime,m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*scale3.*epsU/scale5);
% %         hold on; plot(MITTime,m0*(365.25*24*60*60)*TminTf_IF.^2.*mu.^2.*epsT.*scale4/scale5);
% 
%         [R,P] = corrcoef(movmean_TminTf,movmean_melt);
%         R2_Tstar2(gg,1) = R(1,2)^2*100;
%         P_Tstar2(gg,1) = P(1,2);
%         scale_Tstar(gg,1) = sqrt(scale1);
%         [R,P] = corrcoef(movmean_mu,movmean_melt);
%         R2_mu2(gg,1) = R(1,2)^2*100;
%         P_mu2(gg,1) = P(1,2);
%         scale_mu(gg,1) = sqrt(scale2);
%         [R,P] = corrcoef(movmean_epsU,movmean_melt);
%         R2_epsU(gg,1) = R(1,2)^2*100;
%         P_epsU(gg,1) = P(1,2);
%         scale_epsU(gg,1) = scale4;
%         [R,P] = corrcoef(movmean_epsT,movmean_melt);
%         R2_epsT(gg,1) = R(1,2)^2*100;
%         P_epsT(gg,1) = P(1,2);
%         scale_epsT(gg,1) = scale3;
%         scale_melt(gg,1) = scale5;
% 

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

%    T = table(basins',scale_Tstar,MAPE_Tstar2,RMSE_Tstar2,scale_mu,MAPE_mu2,RMSE_mu2,scale_epsT,MAPE_epsT,RMSE_epsT,scale_epsU,MAPE_epsU,RMSE_epsU,scale_melt)


end


for kk=1:4

    patch(ax_fig(kk),'Xdata',[datenum('01011980','ddmmyyyy') datenum('01011980','ddmmyyyy') ...
        datenum('01012015','ddmmyyyy') datenum('01012015','ddmmyyyy')],'Ydata',[-10 100 100 -10],'facecolor',[0.5 0.5 0.5],'facealpha',0.2,'linestyle','--','linewidth',2,'Edgecolor',[0.6 0.6 0.6]);
    xlim(ax_fig(kk),[datenum('01011990','ddmmyyyy') datenum('01012215','ddmmyyyy')]); 
    xticks(ax_fig(kk),[datenum('01012015','ddmmyyyy') datenum('01012040','ddmmyyyy') datenum('01012065','ddmmyyyy')...
        datenum('01012090','ddmmyyyy') datenum('01012115','ddmmyyyy') datenum(['01012140'],'ddmmyyyy') ...
        datenum('01012165','ddmmyyyy') datenum('01012190','ddmmyyyy') datenum('01012215','ddmmyyyy')]);
    xticklabels(ax_fig(kk),["0","","50","","100","","150","","200"]);

    title(ax_fig(kk),"$\textbf{"+basinslegend(kk)+"}$",'interpreter','latex','fontsize',14);
    %h(1) = plot([0 1],[-1 -1],'-k','linewidth',2);
    %h(2) = plot([0 1],[-1 -1],'--k','linewidth',2);
    
    grid(ax_fig(kk),'on'); box(ax_fig(kk),'on');    

end

xticklabels(ax_fig(1),{});
xticklabels(ax_fig(2),{});
%yyaxis(ax_fig(1),"left");
ylim(ax_fig(1),[20 95]);
%yyaxis(ax_fig(2),"left");
ylim(ax_fig(2),[30 60]);
%yyaxis(ax_fig(3),"left");
ylim(ax_fig(3),[10 45]);
%yyaxis(ax_fig(4),"left");
ylim(ax_fig(4),[5 35]);


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
% 
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
% 
% for jj=1:numel(runID)
%     for gg=[1 3 4 2]
%         uistack(g1((jj-1)*4+gg),'top');
%         uistack(g2((jj-1)*4+gg),'top');
%         uistack(g3((jj-1)*4+gg),'top');
%         uistack(g4((jj-1)*4+gg),'top');
%         uistack(g5((jj-1)*4+gg),'top');
% %        uistack(g6(gg),'top');
% %        uistack(g8(gg),'top');
%     end
% end

%%leg1=legend(g1(1:4),basinslegend{:},'Orientation', 'Horizontal','fontsize',12);
%%leg1.Layout.Tile = 'north';

%ah1=axes('position',get(gca,'position'),'visible','off');

%leg2=legend(ah1,h(1:2),{'$\mbox{\textit{hi\_melt}}$','$\mbox{\textit{av\_melt}}$'},...
%    'Orientation', 'Horizontal','fontsize',14,'interpreter','latex','FontWeight','bold');
%leg2.Position = [0.1 0.95 0.14 0.02];

xlabel(tlo_fig,'$\rm{time \,[yr]}$','interpreter','latex','fontsize',14);
ylabel(tlo_fig,"$\rm{Average\,melt\; [m\, yr}^{-1}\rm{]}$",'interpreter','latex','fontsize',14);

%yyaxis(ax_fig(2),"right")
%ylabel(ax_fig(2),"$\epsilon_U \; \rm{[ms}^{-1\circ}\rm{C}^{-1}\rm{]}$","interpreter","latex",'FontSize',14);
%yyaxis(ax_fig(4),"right")
%ylabel(ax_fig(4),"$\epsilon_U$","interpreter","latex");
% 
%leg=legend(ax_fig(1),[g2(3) g2(1) g1(1) g2(2) g1(2) ],["$m \; \mbox{(\textit{av\_melt})}$","$m \; \mbox{(\textit{ref\_melt})}$","$\epsilon_U \; \mbox{(\textit{ref\_melt})}$",...
%    "$m \; \mbox{(\textit{var\_melt})}$","$\epsilon_{U} \; \mbox{(\textit{var\_melt})}$"],'Orientation', 'Horizontal',"interpreter","latex","fontsize",14);
%leg.Layout.Tile = 'north';
leg=legend(ax_fig(1),[g2(:)],["$\mbox{\textit{ref\_melt}}$","$\mbox{\textit{var\_melt}}$","$\mbox{\textit{av\_melt}}$"],'Orientation', 'Horizontal',"interpreter","latex","fontsize",14);
%leg.Layout.Tile = 'north';
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
fname = ['../Figures/Fig7_GeometricCoefficients_varmelt'];
print(H1,fname,'-dpng','-r400');

for gg=1:numel(basins)
    basin = char(basins(gg));
    diff = mean(meltdata(2).(basin).movmean_melt-meltdata(1).(basin).movmean_melt,"all","omitnan");
    pdiff = 100*diff/mean(meltdata(1).(basin).movmean_melt,"all","omitnan");
    fprintf('====%s====\n',char(basin));
    fprintf('Mean difference ref - var: %s\n',num2str(mean(diff,"all","omitnan")));
    fprintf('Percentage mean difference 100*(ref - var)/ref: %s\n',num2str(pdiff));
end

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