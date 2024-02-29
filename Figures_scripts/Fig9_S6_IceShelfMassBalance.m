function Fig9_S6_IceShelfMassBalance

addpath(getenv("froot_tools"));

%% components
%1. dh/dt
%2. GL flux
%3. Icefront flux
%4. surface accumulation
%5. basal ablation

runID= "PTDC_001";

%years = [2017:2250];

basins=DefineBasins;
basintitle = ["Pine Island","Thwaites","Crosson","Dotson"];
    
if ~exist("IceShelfMassBalance_"+runID+".mat","file")

    frootm = getenv('froot_uamitgcm')+"cases/"+runID;
    subd=dir(frootm+"/output/");
    isub = [subd(:).isdir]; %# returns logical vector
    nameFolds = {subd(isub).name}';
    nameFolds(ismember(nameFolds,{'.','..'})) = [];
    
    output = []; kk=1;
    switch runID
        case "PTDC_001"
            endtime=2726;
        case {"PTDC_002","PTDC_003"}
            endtime = 12*200;
    end

    for jj=1:endtime
        
        output(jj).filelist=dir(frootm+"/output/"+nameFolds{jj}+"/Ua/UaDefault*.mat");
        for ii=1:length(output(jj).filelist)
            files{kk}.folder = output(jj).filelist(ii).folder;
            files{kk}.name = output(jj).filelist(ii).name;
            kk = kk+1;
        end
    end
    
    kk=1;
    while kk<=numel(files)
    
        files{kk}.name
        
        Uafile = files{kk}.folder+"/"+files{kk}.name;
        nstr = strlength(Uafile);
        yyyymm = extractBetween(Uafile,nstr-14,nstr-9);
        dd = double(extractBetween(Uafile,nstr-7,nstr-4));
        data.time(kk) = datenum(yyyymm+string(dd),"yyyymmdd");    

        load(Uafile);
    
        [xB,yB,qB,IceFrontNodes]=FluxAcrossModelBoundary(CtrlVar,MUA,ub,vb,h,b,B,rho); %kg/yr
        xB = 1/2*(xB(:,1)+xB(:,2));
        yB = 1/2*(yB(:,1)+yB(:,2));
        qB_IF = qB; qb_IF(~IceFrontNodes)=0;
        [GF,~,~,~]=IceSheetIceShelves(CtrlVar,MUA,GF);
        %[~,LakeNodes,~,~,~,~] = LakeOrOcean(CtrlVar,MUA,GF);
        %[LakeNodes,OceanNodes,LakeElements,OceanElements] = LakeOrOcean3(CtrlVar,MUA,GF);
        [qGL,qGLx,qGLy,Fub,Fvb,Fr,Fh,LakeNodes,GLgeo,ubGL,vbGL] = FluxAcrossGroundingLine_JDR(CtrlVar,MUA,GF,ub,vb,0*ub,0*vb,h,rho); %kg/yr
        ISArea_int = FEintegrate2D(CtrlVar,MUA,1-GF.node);
        ISArea_below400 = 1-GF.node; ISArea_below400(b>-400)=0;
        ISArea_below400_int = FEintegrate2D(CtrlVar,MUA,ISArea_below400);
        as_int =  FEintegrate2D(CtrlVar,MUA,as*1e3);
        as_ISint = FEintegrate2D(CtrlVar,MUA,1e3*as.*(1-GF.node)); %kg/yr
        %as_ISint(LakeElements)=0;
        ab_int = FEintegrate2D(CtrlVar,MUA,ab*1e3); %kg/yr
        ab_below400 = ab; ab_below400(b>-400)=0;
        ab_below400_int = FEintegrate2D(CtrlVar,MUA,ab_below400*1e3); %kg/yr
        %ab_int(LakeElements)=0;
        dhdt_int  = FEintegrate2D(CtrlVar,MUA,dhdt.*rho);
        dhdt_ISint = FEintegrate2D(CtrlVar,MUA,dhdt.*rho.*(1-GF.node)); %kg/yr
        %dhdt_ISint(LakeElements)=0;
        xEle = Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,1));
        yEle = Nodes2EleMean(MUA.connectivity,MUA.coordinates(:,2));
       
        for bb=1:numel(basins)
            
            basin = basins(bb).Name;
            xpoly = basins(bb).X;
            ypoly = basins(bb).Y;
            
            IIF = find(inpoly([xB(:) yB(:)],[xpoly(:) ypoly(:)]));
            IGL = find(inpoly([GLgeo(:,7),GLgeo(:,8)],[xpoly(:) ypoly(:)]));
            IEle = find(inpoly([xEle(:) yEle(:)],[xpoly(:) ypoly(:)]));
    
            %0. Ice Shelf Area
            data.(basin).A(kk) = sum(ISArea_int(IEle),"all"); %m^2
            data.(basin).A_below400(kk) = sum(ISArea_below400_int(IEle),"all"); %m^2
            %1. dh/dt
            data.(basin).dhdt(kk) = sum(dhdt_int(IEle),"all"); %m^3/yr to kg/yr
            data.(basin).ISdhdt(kk) = sum(dhdt_ISint(IEle),"all");
            %2. GL flux
            data.(basin).GLflux(kk) = sum(qGL(IGL),"all"); %kg/yr
            %3. Icefront flux
            data.(basin).IFflux(kk) = sum(qB(IIF),"all"); %kg/yr
            %4. surface accumulation
            data.(basin).as(kk) = sum(as_int(IEle),"all");
            data.(basin).ISas(kk) = sum(as_ISint(IEle),"all"); %kg/yr
            %5. basal ablation
            data.(basin).ab(kk) = sum(ab_int(IEle),"all"); %kg/yr
            data.(basin).ab_below400(kk) = sum(ab_below400_int(IEle),"all"); %kg/yr
            %6. sum
            data.(basin).sum(kk) = data.(basin).ab(kk)+data.(basin).as(kk)+data.(basin).IFflux(kk)-data.(basin).dhdt(kk); %kg/yr
            data.(basin).ISsum(kk) = data.(basin).ab(kk)+data.(basin).ISas(kk)+data.(basin).IFflux(kk)+data.(basin).GLflux(kk)-data.(basin).ISdhdt(kk); %kg/yr
    
        end
    
        kk=kk+1;
    
    end

    save("IceShelfMassBalance_"+runID+".mat","data");

else

    load("IceShelfMassBalance_"+runID+".mat");

end

%% PLOTTING
H=fig('units','inches','width',90*12/72.27,'height',45*12/72.27,'fontsize',14,'font','Helvetica');

tlo_fig = tiledlayout(1,4,"TileSpacing","tight");
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

kk=1;


switch runID
    case {"PTDC_002","PTDC_003"}
        data2=load("IceShelfMassBalance_PTDC_001.mat");
        for gg=1:numel(basins)
            basin = basins(gg).Name;
            data.(basin).A = [data2.data.(basin).A(1:18*12) data.(basin).A];
            data.(basin).ab = [data2.data.(basin).ab(1:18*12) data.(basin).ab];
            data.(basin).ab_below400 = [data2.data.(basin).ab_below400(1:18*12) data.(basin).ab_below400];
            data.(basin).GRas = [data2.data.(basin).as(1:18*12)-data2.data.(basin).ISas(1:18*12) data.(basin).as-data.(basin).ISas];
            data.(basin).ISas = [data2.data.(basin).ISas(1:18*12) data.(basin).ISas];
            data.(basin).IFflux = [data2.data.(basin).IFflux(1:18*12) data.(basin).IFflux];
            [~,IqGLmax(gg)]=max(data.(basin).GLflux);
            data.(basin).GLflux = [data2.data.(basin).GLflux(1:18*12) data.(basin).GLflux];
            timeGLmax(gg) = (data.time(IqGLmax(gg))-data.time(1))/365.25;
        end
        time = ([data2.data.time(1:18*12) data.time]-data.time(1))/365.25;
        
    case "PTDC_001"
        for gg=1:numel(basins)
            basin = basins(gg).Name;
            [~,IqGLmax(gg)]=max(data.(basin).GLflux);
            timeGLmax(gg) = (data.time(IqGLmax(gg))-data.time(18*12))/365.25;
            data.(basin).GRas = data.(basin).as-data.(basin).ISas;
        end
        time = (data.time-data.time(18*12))/365.25;
end

for gg=2:5

    basin=basins(gg).Name;

    g(7)=patch(ax_fig(kk),'Xdata',[-30 -30 0 0],'Ydata',[-200 200 200 -200],'facecolor',[0.8 0.8 0.8],'facealpha',0.6,'linestyle','none');
    
%     if timeGLmax(gg)>5
%         plot(ax_fig(kk),[timeGLmax(gg) timeGLmax(gg)],[-200 200],'--','linewidth',2,'color',[0.65 0.65 0.65]);
%     end

    %yyaxis(ax_fig(kk),"left");
    plot(ax_fig(kk),time,data.(basin).ab_below400/1e12,'-','color',[0.5 0.9 0.97]);
    
    
    plot(ax_fig(kk),time,(data.(basin).ab-data.(basin).ab_below400)/1e12,'--','color',[0.5 0.9 0.97]);
    dh=data.(basin).GLflux+data.(basin).ISas+data.(basin).ab+data.(basin).IFflux;
    plot(ax_fig(kk),time,dh/1e12,'-','color',[0.7 0.7 0.7]);%-dh(1)); 

    g(8)=plot(ax_fig(kk),time,data.(basin).GRas/1e12,'--','color',[0.47 0.67 0.19],'linewidth',2);%-data.(basin).ISas(1)); 
    g(4)=plot(ax_fig(kk),time,data.(basin).ISas/1e12,'-','color',[0.47 0.67 0.19],'linewidth',2);%-data.(basin).ISas(1)); 
    g(5)=plot(ax_fig(kk),time,data.(basin).IFflux/1e12,'-','color',[0.93 0.69 0.13],'linewidth',2);%-data.(basin).IFflux(1));
    g(6)=plot(ax_fig(kk),time,data.(basin).GLflux/1e12,'-','color',[0.64 0.08 0.18],'linewidth',2);%-data.(basin).GLflux(1));   

    g(1)=plot(ax_fig(kk),time,movmean(data.(basin).ab_below400/1e12,5*12),'-','color',[0 0.45 0.74],'linewidth',2);
    g(2)=plot(ax_fig(kk),time,movmean((data.(basin).ab-data.(basin).ab_below400)/1e12,5*12),'--','color',[0 0.45 0.74],'linewidth',2);
    g(3)=plot(ax_fig(kk),time,movmean(dh/1e12,5*12),'-','color',[0 0 0],'linewidth',2);
    %patch(ax_fig(kk),'Xdata',[-30 -30 0 0],'Ydata',[-200 200 200 -200],'facecolor',[1 1 1],'facealpha',0.6,'linestyle','none');
    %plot(ax_fig(kk),[timeGLmax(gg) timeGLmax(gg)],[-200 200],'--k','linewidth',2)
    
    
%     yyaxis(ax_fig(kk),"right");
%     plot(ax_fig(kk),time,movmean(data.(basin).ab/1e3./data.(basin).A,12*5),'-k');
%     ylim(ax_fig(kk),[-20 0])

    %yyaxis(ax_fig(kk),"left");
    xticks(ax_fig(kk),[0 50 100 150 200]);
    xticklabels(ax_fig(kk),["0" "" "100" "" "200"]);
    xlim(ax_fig(kk),[-20 200]);

    switch kk
        case 1
            ylim(ax_fig(kk),[-120 120]);
            ylabel(tlo_fig,'$\mbox{Mass balance components [Gt yr}^{-1}\mbox{]}$','interpreter','latex','fontsize',14);
            xlabel(tlo_fig,'$\mbox{Time [yr]}$','interpreter','latex','fontsize',14);
        case 2
            ylim(ax_fig(kk),[-120 120]);
        case 3
            ylim(ax_fig(kk),[-45 45]);
        case 4
            ylim(ax_fig(kk),[-45 45]);
            %lgd=legend([g(6) g(5) g(3) g(1) g(4) g(2)],["Spinup+Hindcast","Grounding line flux","Surface mass balance","Basal melt","Calving flux","Thickness change"]);
            lgd=legend([g(7) g(6) g(8) g(4) g(5) g(1) g(2)  g(3)],["Hindcast 1997-2014","Grounding line flux","SMB (grounded ice)","SMB (ice shelf)","Calving flux","Ice-shelf melt (deep interior)","Ice-shelf melt (outer cavity)","Ice-shelf mass change"]);
            lgd.Orientation="horizontal";
            lgd.Layout.Tile ="north";
            lgd.NumColumns=4;
    end
    grid(ax_fig(kk),"on"); box(ax_fig(kk),"on"); 
    title(ax_fig(kk),basintitle(kk))

    kk=kk+1;

end

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
folder = ["../Figures/"];
switch runID
    case "PTDC_002"
        fname = [folder+"/Fig8_IceShelfMassBalance_"+runID];
    case "PTDC_001"
        fname = [folder+"/FigS6_IceShelfMassBalance_"+runID];
end
print(H,fname,"-dpng","-r400");


end

%%%%%%%%%%%%%%%%%%%
%% BASIN OUTLINE %%
%%%%%%%%%%%%%%%%%%%
function [xpoly,ypoly] = BasinOutline(basin)

    switch basin
        case "PIG"
            xpoly = [-1699 -1699 -1500 -1500 -1550]*1e3;
            ypoly = [-375 -220 -220 -300 -375]*1e3;
        case "TW"
            xpoly = [-1620 -1550 -1500 -1450 -1450 -1620]*1e3;
            ypoly = [-375 -375 -300 -300 -520 -520]*1e3;
        case "CR"
            xpoly = [-1610 -1485 -1450 -1450 -1610]*1e3;
            ypoly = [-580 -657 -657 -520 -520]*1e3;
        case "DT"
            xpoly = [-1610 -1485 -1450 -1450 -1625]*1e3;
            ypoly = [-580 -657 -657 -700 -700]*1e3; 
        case "AS"
            xpoly = [-2e6 -2e6 -0.5e6 -0.5e6]; 
            ypoly = [-9e5 1.1e5 1.1e5 -9e5];
    end

end

    