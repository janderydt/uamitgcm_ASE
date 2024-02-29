function Fig7_StreamFunctionAmplitude

runID = "PTDC_002";
basins = ["PIG","TW","CR","DT"];
glaciername = ["Pine Island","Thwaites","Crosson","Dotson"];

%color = [[51 34 136];[17 119 51];[133 87 35];[136 34 85]]/255;
%color = [[31,120,180];[27,158,119];[246,190,0];[152,78,163]]/255;
%colorsoft = [[166,206,227];[141,211,199];[251,219,101];[190,174,212]]/255;

%color = [31,120,180]/255;
%colorsoft = [166,206,227]/255;
%colorepsU = [178,223,138]/255;
btvelcolor = [0.47 0.67 0.19];
bsfcolor = [0.75 0.33 0.1];
bsfcolor2 = [0.89 0.74 0.44];
osfcolor = [0.91 0.69 0.13];

folder = "../Figures/"+runID;
if ~exist(folder,"dir")
    mkdir(folder);
end

froot = getenv("froot_uamitgcm");

frootm = froot+"cases/"+runID;

load(froot+"/Ua_InputData/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat","FB");

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = string({subd(isub).name});
nameFolds(ismember(nameFolds,[".",".."])) = [];
	
output = []; kk=1;
for jj=1:numel(nameFolds)
    disp(frootm+"/output/"+nameFolds(jj));
    output(jj).filelist=dir(frootm+"/output/"+nameFolds(jj)+"/MITgcm/output.nc");
    for ii=1:length(output(jj).filelist)
        listing{kk}.folder = output(jj).filelist(ii).folder;
        listing{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
end

start = 1; kk=1;

if exist(getenv("froot_tools")+"CavityStreamFunctions_"+runID+"_monthly_v2.mat","file")~=2

    for ii=start:start+200*12

        addpath(getenv("froot_tools"));

	    MITfile=listing{ii}.folder+"/"+listing{ii}.name;

        nstr = strlength(MITfile);
        yyyymm = extractBetween(MITfile,nstr-22,nstr-17);
        time(kk) = datenum(yyyymm+"01","yyyymmdd");        
    
        [bsf_tmp,osf_tmp] = CalcCavityStreamFunctions(runID,MITfile,1);
        
    
        for gg=1:numel(basins)

            basin = basins{gg};

            bsf.(basin).fixed(kk).amp = max(bsf_tmp.(basin).fixedmask(:))-min(bsf_tmp.(basin).fixedmask(:));
            bsf.(basin).fixed(kk).max = max(bsf_tmp.(basin).fixedmask(:));
            bsf.(basin).moving(kk).amp = max(bsf_tmp.(basin).movingbelow400m(:))-min(bsf_tmp.(basin).movingbelow400m(:));
            bsf.(basin).moving(kk).max = max(bsf_tmp.(basin).movingbelow400m(:));
            bsf.(basin).moving(kk).mask = bsf_tmp.maskbelow400m.(basin);

	        osf.(basin).amp(kk) = max(osf_tmp.(basin)(:))-min(osf_tmp.(basin)(:));
	        osf.(basin).max(kk) = max(osf_tmp.(basin)(:));
	    
        end

	    fprintf("Done %i out of %i\n",ii,200*12);	
	    kk = kk+1;	
	    
    end

    for gg=1:numel(basins)

        basin = basins{gg};
        bsf.(basin).fixedmask = bsf_tmp.fixedmask.(basin);

    end
    
    save(getenv("froot_tools")+"CavityStreamFunctions_"+runID+"_monthly_v2.mat","time","bsf","osf","-v7.3");

else

    load(getenv("froot_tools")+"CavityStreamFunctions_"+runID+"_monthly_v2.mat");
    SFtime = time;
    load(getenv("froot_tools")+"HeatVolumeTransport_moving400mdraft_"+runID+".mat");

     for ii=[start+2*12 200*12]

        addpath(getenv("froot_tools"));

	    MITfile=listing{ii}.folder+"/"+listing{ii}.name;
        nstr = strlength(MITfile);
        folder = extractBetween(MITfile,1,nstr-17);
        Uafile = dir(folder+"/Ua/UaDefault*.mat");
        load(Uafile.folder+"/"+Uafile.name);
        CtrlVar.PlotGLs=0;
       
        nstr = strlength(MITfile);
        yyyymm = extractBetween(MITfile,nstr-22,nstr-17);
        time = datenum(yyyymm+"01","yyyymmdd");        

        [xGL.(char("bsf_"+yyyymm+"01")),yGL.(char("bsf_"+yyyymm+"01")),GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF);
        xB.(char("bsf_"+yyyymm+"01")) = MUA.Boundary.x;
        yB.(char("bsf_"+yyyymm+"01")) = MUA.Boundary.y;

        [bsf_tmp,osf_tmp] = CalcCavityStreamFunctions(runID,MITfile,1);

        bsf.AS.(char("bsf_"+yyyymm+"01")) = bsf_tmp.AS.fixedmask;
        
        Lon = double(ncread(MITfile,'XC'));
        Lat = double(ncread(MITfile,'YC'));
        if size(Lon,2)>1
            Lon = Lon(:,1);
            Lat = Lat(1,:);
        end
        [X,Y] = ndgrid(Lon,Lat);
        bsf.AS.(char("X_"+yyyymm+"01")) = X;
        bsf.AS.(char("Y_"+yyyymm+"01")) = Y;

        [~,I] = min(abs(time-MITTime));
        bsf.(char("contours_"+yyyymm+"01")).xmid = section(I).xmid;
        bsf.(char("contours_"+yyyymm+"01")).ymid = section(I).ymid;

     end

     endtime = datenum('01012215','ddmmyyyy');
     [~,end_MITTime] = min(abs(MITTime - endtime));
     kk=1;
     for gg=1:numel(basins)

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
        end
        TminTf=[]; Ustar_bl=[];
         for tt=1:end_MITTime
    
                I = find(inpoly([section(tt).xmid(:) section(tt).ymid(:)],[xpoly(:) ypoly(:)]));
                
                % thermal driving, heat flux, volume flux, cross section
                %TminTf(kk) = mean(section(tt).monthly.TminTf_belowz_in_depthmean(I),'all','omitnan');
                TminTf(tt) = mean(section(tt).monthly.TminTf_in_depthmean(I),'all','omitnan');
                Ustar_bl(tt) = integral2D.(basin).monthly.blUStar_mean(tt);
                melt_tmp(tt) = integral2D.(basin).monthly.melt_integral(tt)./integral2D.(basin).monthly.ISarea_integral(tt)*1e12/1e3; %m/yr;

         end

         epsU_tmp  = Ustar_bl./TminTf;
         epsU.(basin) = interp1(MITTime(1:end_MITTime),epsU_tmp,SFtime);
         melt.(basin) = interp1(MITTime(1:end_MITTime),melt_tmp,SFtime);

     end

end

H=fig('units','inches','height',60*12/72.27,'width',120*12/72.27,'fontsize',14,'font','Helvetica','visible','off');

%set(H, "visible","off");
%tlo = tiledlayout(2,4,"TileSpacing","tight");

time = (SFtime-SFtime(1))/365.25;

%nexttile([2,1]);
subplot("position",[0.59 0 0.41 1]); hold on;

MITfile=listing{start}.folder+"/"+listing{start}.name;
Melt = double(ncread(MITfile,"SHIfwFlx"));
Mask = Melt; Mask(Mask~=0)=1; Mask(Mask==0)=NaN;
Lon = double(ncread(MITfile,"XC"));
Lat = double(ncread(MITfile,"YC"));
if size(Lon,2)>1
    Lon = Lon(:,1);
    Lat = Lat(1,:);
end
[X,Y]=ndgrid(Lon,Lat); 
%pcolor(X,Y,bsf.AS.bsf_20170101.*Mask); shading flat;
contourf(X,Y,bsf.AS.bsf_20170101.*Mask,[-2:0.5:8]*1e5,'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5); 
CM1=othercolor('RdYlBu_11b',4*32+30); CM1 = flipdim(CM1,1);
CM2=othercolor('RdYlBu_11b',16*32); CM2 = flipdim(CM2,1);
CM=[CM1(1:65,:); CM2(7*32:end,:)];
colormap(CM);
cb=colorbar("Position",[0.9246,0.5416,0.007,0.35],"Location","east","AxisLocation","out","Ticks",[-0.2:0.2:0.8]*1e6,...
    "TickLabels",["-0.2" "0" "0.2" "0.4" "0.6" "0.8"]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.04;
cb.FontSize=14;
cb.Label.String = "Barotropic Streamfunction [Sv]";
caxis([-2 8]*1e5);
hold on;
plot(xB.bsf_20170101,yB.bsf_20170101,'-k');
plot(xGL.bsf_20170101,yGL.bsf_20170101,'-','color',[0.2 0.2 0.2],'linewidth',1);
plot(bsf.contours_20170101.xmid,bsf.contours_20170101.ymid,'.m','markersize',4);

% Lon = Lon(:); Lat = Lat(:);
% %PIG
% plot([Lon(58) Lon(58) Lon(110) Lon(110) Lon(58)]+0.18e6,[Lat(310) Lat(end-20) Lat(end-20) Lat(310) Lat(310)],":","color",colorsoft,"linewidth",3);
% %TW
% plot([Lon(110) Lon(110) Lon(147) Lon(147) Lon(110)]+0.18e6,[Lat(150) Lat(260) Lat(260) Lat(150) Lat(150)],":","color",colorsoft,"linewidth",3);
% %CR
% plot([Lon(130) Lon(130) Lon(190) Lon(190) Lon(130)]+0.18e6,[Lat(35) Lat(110) Lat(110) Lat(35) Lat(35)],":","color",colorsoft,"linewidth",3);
% %DT
% plot([Lon(140) Lon(140) Lon(190) Lon(190) Lon(140)]+0.18e6,[Lat(5) Lat(33) Lat(33) Lat(5) Lat(5)],":","color",colorsoft,"linewidth",3);

axis equal; 
ylim([-7e5 -2.25e5]); xlim([-1.69e6 -1.43e6]);
xticklabels(gca,{}); yticklabels(gca,{});
set(gca,"Visible","off");

%% depth average flow
% quiver(X,Y,u_bt.*xymask,v_bt.*xymask,"color","k","AutoScaleFactor",10);

%subplot("position",[0.251 0 0.4 1]); hold on;

MITfile=listing{start+12*200}.folder+"/"+listing{start+12*200}.name;
Melt = double(ncread(MITfile,"SHIfwFlx"));
Mask = Melt; Mask(Mask~=0)=1; Mask(Mask==0)=NaN;
Lon = double(ncread(MITfile,"XC"));
Lat = double(ncread(MITfile,"YC"));
if size(Lon,2)>1
    Lon = Lon(:,1);
    Lat = Lat(1,:);
end
[X,Y]=ndgrid(Lon,Lat); 
contourf(X+0.16e6,Y,bsf.AS.bsf_22141201.*Mask,[-2:0.5:8]*1e5,'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.5); 
caxis([-2 8]*1e5);
hold on
I = find(xB.bsf_20170101<-1.6894e6); 
xB.bsf_20170101(I)=NaN; yB.bsf_20170101(I)=NaN;
plot(xB.bsf_20170101+0.16e6,yB.bsf_20170101,'-k');
%plot(xGL.bsf_20170101,yGL.bsf_20170101,'-k');
plot(xGL.bsf_22141201+0.16e6,yGL.bsf_22141201,'-','color',[0.2 0.2 0.2],'linewidth',1);
plot(bsf.contours_22141201.xmid+0.16e6,bsf.contours_22141201.ymid,'.m','markersize',4);

Lon = Lon(:); Lat = Lat(:);
%PIG
plot([Lon(58) Lon(58) Lon(110) Lon(110) Lon(58)]+0.16e6,[Lat(310) Lat(end-20) Lat(end-20) Lat(310) Lat(310)],":","color",bsfcolor2,"linewidth",3);
%TW
plot([Lon(110) Lon(110) Lon(147) Lon(147) Lon(110)]+0.16e6,[Lat(150) Lat(260) Lat(260) Lat(150) Lat(150)],":","color",bsfcolor2,"linewidth",3);
%CR
plot([Lon(130) Lon(130) Lon(190) Lon(190) Lon(130)]+0.16e6,[Lat(35) Lat(110) Lat(110) Lat(35) Lat(35)],":","color",bsfcolor2,"linewidth",3);
%DT
plot([Lon(140) Lon(140) Lon(190) Lon(190) Lon(140)]+0.16e6,[Lat(5) Lat(33) Lat(33) Lat(5) Lat(5)],":","color",bsfcolor2,"linewidth",3);

axis equal; ylim([-7e5 -2.25e5]); xlim([-1.69e6 -1.45e6+0.2e6]);
xticklabels(gca,{}); yticklabels(gca,{});
set(gca,"Visible","off");

for gg=1:numel(basins)

    basin = basins{gg};

    ax1(gg) = axes("position",[0.05+(1-mod(gg,2))*0.235 0.1+(1-floor(gg/3))*0.455 0.218 0.38]); hold on;

    yyaxis(ax1(gg),'left');
    %plot(ax1(gg),time,epsU.(basin),"-","color","k","linewidth",1.5);
    plot(ax1(gg),time,melt.(basin),"-","color","k","linewidth",1.5);
    xlim(ax1(gg),[0 200]);
    
    yyaxis(ax1(gg),'right');
    I = find([bsf.(basin).moving(:).amp]~=0);
    h = plot(ax1(gg),time(I),integral2D.(basin).monthly.btVEL_mean(I),"-","color",btvelcolor,"linewidth",1.5);

    ax2(gg) = axes('position', ax1(gg).Position); hold on;

    yyaxis(ax2(gg),'left');
    ax2(gg).Color = 'none';
    grid(ax2(gg), 'off');
    ax2(gg).XLim = ax1(gg).XLim;
    ax2(gg).XTick = ax1(gg).XTick;

    yyaxis(ax2(gg),'right');
    %ax2(gg).Color = 'none';
    %grid(ax2(gg), 'off');
    %ax2(gg).YLimMode = 'manual';
    g(1)=plot(ax2(gg),time(I),[bsf.(basin).moving(I).amp]/1e6,"-","color",bsfcolor,"linewidth",1);
    g(2)=plot(ax2(gg),time(I),[osf.(basin).amp(I)]/1e6,"-","color",bsfcolor2,"linewidth",1);

    tmp = epsU.(basin); 
    %[R,P] = corrcoef(tmp(I),[bsf.(basin).moving(I).amp]/1e6,'rows','complete','alpha',0.05);
    [R,P] = corrcoef(tmp(I),integral2D.(basin).monthly.btVEL_mean(I),'rows','complete');
    fprintf('Basin: %s - Correlation between eps_U and average depth-mean flow: %d \n',basin,R(1,2));
    
    [R2,P2] = corrcoef([bsf.(basin).moving(I).amp],[osf.(basin).amp(I)],'rows','complete');
    fprintf('Basin: %s - Correlation between barotropic & baroclinic streamfunction amplitude: %d \n',basin,R2(1,2));

    grid(ax1(gg),'on');
    box(ax1(gg),'on');

    ax1(gg).YAxis(1).Color = [0 0 0];
    ax1(gg).YAxis(2).Color = btvelcolor;
    ax2(gg).YAxis(1).Color = 'none';
    ax2(gg).XAxis.Color = 'none';
    ax2(gg).YAxis(2).Color = 'none';

    
    if gg==1
        xticklabels(ax1(gg),{});
        
        yyaxis(ax1(gg),'left');
        %ylabel(ax1(gg),"$\epsilon_U \; \rm{[ms}{}^{-1 \circ}\rm{C}{}^{-1}\rm{]}$","interpreter","latex","FontSize",16);        
        ylabel(ax1(gg),"melt rate [m/yr]");
        %ylim(ax1(gg),[0 3.7]*1e-3);
        %yticks(ax1(gg),[0:0.5:3.5]*1e-3);
        ylim(ax1(gg),[0 120]);
        yticks(ax1(gg),[0:20:120]);

        text(ax1(gg),0.95,0.05,['R = ',num2str(R(1,2),2)],Units='normalized',HorizontalAlignment='right',VerticalAlignment='bottom',Margin=5,EdgeColor='k',BackgroundColor='w');

        yyaxis(ax2(gg),'left');
        ax2(gg).YAxis(1).Color='none';
        yticklabels(ax2(gg),{});

        yyaxis(ax1(gg),'right');
        yticklabels(ax1(gg),{});
        ylim(ax1(gg),[0 2]*1e-1);
        yticks(ax1(gg),[0:0.25:2]*1e-1);
        
        yyaxis(ax2(gg),'right');
        yticklabels(ax2(gg),{});
        ylim(ax2(gg),[0 12]*1e5/1e6);
        
    elseif gg==2
        xticklabels(ax1(gg),{});

        ax3(gg) = axes('position', ax1(gg).Position+[0.06 0 0 0]); hold on;
        ax3(gg).XAxis.Color='none';
        yyaxis(ax3(gg),'right');
        ax3(gg).Color = 'none';
        grid(ax3(gg), 'off');
        ax3(gg).XLim = ax1(gg).XLim;
        ax3(gg).XTick = ax1(gg).XTick;
        ax3(gg).YAxis(1).Color='none';
        ax3(gg).YAxis(2).Color=bsfcolor;
        grid(ax3(gg),'off');

        yyaxis(ax1(gg),'left');
        %ylabel("$\epsilon_U \; \rm{[ms]}{}^{-1 \circ}\rm{C}{}^{-1}\rm{]}$","interpreter","latex");
        yticklabels(ax1(gg),{});
        %ylim(ax1(gg),[0 3.7]*1e-3);
        %yticks(ax1(gg),[0:0.5:3.5]*1e-3);
        ylim(ax1(gg),[0 120]);
        yticks(ax1(gg),[0:20:120]);

        text(ax1(gg),0.05,0.95,['R = ',num2str(R(1,2),2)],Units='normalized',HorizontalAlignment='left',VerticalAlignment='top',Margin=5,EdgeColor='k',BackgroundColor='w');

        yyaxis(ax2(gg),'left');
        yticklabels(ax2(gg),{});

        yyaxis(ax1(gg),'right');
        ylim(ax1(gg),[0 2]*1e-1);
        yticks(ax1(gg),[0:0.25:2]*1e-1);
        
        yyaxis(ax2(gg),'right');
        ylim(ax2(gg),[0 12]*1e5/1e6);
        yticklabels(ax2(gg),{});

        yyaxis(ax3(gg),'right');
        ylim(ax3(gg),[0 12]*1e5/1e6);

        legend(ax1(gg),g(:),{'Barotropic streamfunction','Overturning streamfunction'},'location','northeast');

    elseif gg==3

        yyaxis(ax1(gg),'left');
        %
        %ylabel("$\epsilon_U \; \rm{[ms]}{}^{-1 \circ}\rm{C}{}^{-1}\rm{]}$","interpreter","latex");
        %ylim(ax1(gg),[0 1.8]*1e-3);
        %yticks(ax1(gg),[0:0.25:1.75]*1e-3);
        ylim(ax1(gg),[0 50]);
        yticks(ax1(gg),[0:10:50]);

        yyaxis(ax2(gg),'left');
        yticklabels(ax2(gg),{});
        ax2(gg).YAxis(1).Color='none';

        yyaxis(ax1(gg),'right');
        yticklabels(ax1(gg),{});
        ylim(ax1(gg),[0 0.75]*1e-1);
        yticks(ax1(gg),[0:0.15:0.75]*1e-1);
        
        yyaxis(ax2(gg),'right');
        yticklabels(ax2(gg),{});
        ylim(ax2(gg),[0 4.5]*1e5/1e6);

        xlabel(ax1(gg),"time [a]");
        text(ax1(gg),0.05,0.05,['R = ',num2str(R(1,2),2)],Units='normalized',HorizontalAlignment='left',VerticalAlignment='bottom',Margin=5,EdgeColor='k',BackgroundColor='w');

        
    elseif gg==4

        ax3(gg) = axes('position', ax1(gg).Position+[0.06 0 0 0]); hold on;
        ax3(gg).XAxis.Color='none';
        yyaxis(ax3(gg),'right');
        ax3(gg).Color = 'none';
        grid(ax3(gg), 'off');
        ax3(gg).XLim = ax1(gg).XLim;
        ax3(gg).XTick = ax1(gg).XTick;
        ax3(gg).YAxis(1).Color='none';
        ax3(gg).YAxis(2).Color=bsfcolor;
        grid(ax3(gg),'off');
        ylabel(ax3(gg),"Streamfunction amplitude [Sv]");

        yyaxis(ax1(gg),'left');
        yticklabels(ax1(gg),{});
        %ylabel(ax1(gg),"$\epsilon_U \; \rm{[ms]}{}^{-1 \circ}\rm{C}{}^{-1}\rm{]}$","interpreter","latex","FontSize",16);
        %ylabel("$\epsilon_U \; \rm{[ms]}{}^{-1 \circ}\rm{C}{}^{-1}\rm{]}$","interpreter","latex");
        %ylim(ax1(gg),[0 1.8]*1e-3);
        %yticks(ax1(gg),[0:0.25:1.75]*1e-3);
        ylim(ax1(gg),[0 50]);
        yticks(ax1(gg),[0:10:50]);

        yyaxis(ax2(gg),'left');
        yticklabels(ax2(gg),{});
        ax2(gg).YAxis(1).Color='none';

        yyaxis(ax1(gg),'right');
        ylim(ax1(gg),[0 0.75]*1e-1);
        yticks(ax1(gg),[0:0.15:0.75]*1e-1);
        ylabel(ax1(gg),"Mean depth-averaged velocity [m/s]");
        
        yyaxis(ax2(gg),'right');
        ylim(ax2(gg),[0 4.5]*1e5/1e6);

        text(ax1(gg),0.95,0.95,['R = ',num2str(R(1,2),2)],Units='normalized',HorizontalAlignment='right',VerticalAlignment='top',Margin=5,EdgeColor='k',BackgroundColor='w');
        
    end

    title(ax1(gg),glaciername(gg));

end

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = ["../Figures/Fig7_StreamFunctions_"+runID];   
print(H,fname,"-dpng","-r400");

return

H=fig("units","inches","height",40*12/72.27,"width",45*12/72.27,"fontsize",14,"font","Helvetica","visible","off");

subplot("position",[0.12 0.1 0.8 0.8]); hold on; 

plot([0 2],[0 2],"-k","color",[0.5 0.5 0.5]);

for gg = 1:numel(basins)

    basin = basins{gg};

    nt = numel([bsf.(basin).fixed(:).max]);
    facecolor = [linspace(1,color(gg,1),ceil(1.3*nt))',...
            linspace(1,color(gg,2),ceil(1.3*nt))',linspace(1,color(gg,3),ceil(1.3*nt))'];
    facecolor = facecolor([end-nt+1:end],:);
    g(gg) = scatter([bsf.(basin).moving(:).amp]/1e6,osf.(basin).amp/1e6,13,facecolor,"filled");

end

xlim([0 1.4]); ylim([0 0.4]);
xlabel("bsf [Sv]"); ylabel("osf [Sv]");

grid on; box on;

legend(g,{"Pine Island","Thwaites","Crosson","Dotson"},"Location","southeast","fontsize",14);

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/StreamFunctionRegression_"+runID;   
print(H,fname,"-dpng","-r400");

load(getenv("froot_tools")+"HeatVolumeTransport_moving400mdraft_"+runID+".mat");

for gg = 1:numel(basins)
    
    basin = basins{gg};

    melt_tmp = [integral2D(:).(basin).monthly.melt_integral]./[integral2D(:).(basin).monthly.ISarea_integral]*1e9; % convert from Gt/yr/m2 to m/yr
    melt.(basin) = interp1(MITTime,melt_tmp,SFtime);

end


H=fig('units','inches','height',70*12/72.27,'width',120*12/72.27,'fontsize',14,'font','Helvetica','visible','off');
%set(H, "visible","off");
tlo = tiledlayout(1,2,"TileSpacing","tight");

nexttile; hold on;

for gg = 1:numel(basins)
    
    basin = basins{gg};

    nt = numel(SFtime);
    facecolor = [linspace(1,color(gg,1),ceil(1.3*nt))',...
            linspace(1,color(gg,2),ceil(1.3*nt))',linspace(1,color(gg,3),ceil(1.3*nt))'];
    facecolor = facecolor([end-nt+1:end],:);
    g(gg) = scatter(melt.(basin),[bsf.(basin).moving(:).amp]/1e6,13,facecolor,"filled");

end

xlabel("Average melt [m/yr]"); ylabel("bsf [Sv]");
legend(g,{"Pine Island","Thwaites","Crosson","Dotson"},"Location","northwest","fontsize",14);

grid on; box on;

nexttile; hold on;

for gg = 1:numel(basins)
    
    basin = basins{gg};

    nt = numel(SFtime);
    facecolor = [linspace(1,color(gg,1),ceil(1.3*nt))',...
            linspace(1,color(gg,2),ceil(1.3*nt))',linspace(1,color(gg,3),ceil(1.3*nt))'];
    facecolor = facecolor([end-nt+1:end],:);
    g(gg) = scatter(melt.(basin),osf.(basin).amp/1e6,13,facecolor,"filled");

end

xlabel("Average melt [m/yr]"); ylabel("osf [Sv]");

grid on; box on;

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/Fig7_StreamFunctions";   
print(H,fname,"-dpng","-r400");


