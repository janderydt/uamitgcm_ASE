function FigS6_DepthAvFlow_PV

runID = "PTDC_002";

Itime = [2*12 20*12 40*12 60*12 120*12 180*12];
year=["0","20","40","60","120","180"];

domain = "PIG"; % PIG, TW, DC, AS

addpath(getenv("/mnt/SSD1/Documents/Projects/MISOMIP2/ProcessData"));
frootm = getenv("froot_uamitgcm")+"/cases/"+runID;
addpath(getenv("froot_tools"));

switch domain
    case "PIG"
        %% PIG
        xmin = -1680; xmax = -1550; 
        ymin = -340; ymax = -240;
        labels(1).x=-1660; labels(1).y=-230; labels(1).name = "Pine Island";
    case "TW"
        %% Thwaites
        xmin = -1620; xmax = -1450;
        ymin = -520; ymax = -380;
        labels(1).x=-1530; labels(1).y=-400; labels(1).name = "Thwaites";
    case "DC"
        %% Dotson Crosson
        xmin = -1650; xmax = -1420;
        ymin = -700; ymax = -500;
        labels(1).x=-1485; labels(1).y=-580; labels(1).name = "Pope";
        labels(2).x=-1482; labels(2).y=-620; labels(2).name = "Smith";
        labels(3).x=-1490; labels(3).y=-685; labels(3).name = "Kohler";
    case "AS"
        %% Amundsen
        xmin = -1699; xmax = -1440;
        ymin = -700; ymax = -220;
        labels(1).x=-1660; labels(1).y=-230; labels(1).name = "";
        labels(2).x=-1520; labels(2).y=-420; labels(2).name = "";
        labels(3).x=-1485; labels(3).y=-560; labels(3).name = "";
        labels(4).x=-1478; labels(4).y=-632; labels(4).name = "";
        labels(5).x=-1485; labels(5).y=-685; labels(5).name = "";
end

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = string({subd(isub).name});
nameFolds(ismember(nameFolds,[".",".."])) = [];
	
% read input 
output = []; kk=1;
for jj=1:min(length(nameFolds),Itime(end))
    disp(frootm+"/output/"+nameFolds(jj));
    output(jj).filelist=dir(frootm+"/output/"+nameFolds(jj)+"/Ua/UaDefault*.mat");
    for ii=1:length(output(jj).filelist)
        listing{kk}.folder = output(jj).filelist(ii).folder;
        listing{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
end

load(getenv("froot_tools")+"HeatVolumeTransport_moving400mdraft_"+runID+".mat");

H=fig("units","inches","height",60*12/72.27,"width",110*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(2,3,'TileSpacing','compact');
for i = 1:6
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end


%% Grab data and plot
for ii=1:numel(Itime)

    Uafile = listing{Itime(ii)}.folder+"/"+listing{Itime(ii)}.name;
    MITpath = listing{Itime(ii)}.folder(1:end-3)+"/MITgcm";
    MITfile = MITpath+"/output.nc"; 
    
    Time = double(ncread(MITfile,"time")); Tii=1;
    attvalue=ncreadatt(MITfile,"time","units");
    if strfind(attvalue,"seconds")
        epoch = erase(attvalue,"seconds since ");
        epochnum = datenum(epoch);
        data(ii).Time = double(epochnum+Time(Tii)/(24*60*60));
    elseif strfind(attvalue,"days")
        epoch = erase(attvalue,"days since ");
        epochnum = datenum(epoch);
        data(ii).Time = double(epochnum+Time(Tii));
    else
        error("I do not recognise the time format in output.nc");
    end

    % PV / depth-average flow
    [data(ii).ZLon,data(ii).ZLat,data(ii).fH,data(ii).ZH,data(ii).PV,data(ii).Lon,data(ii).Lat,~,~,~,~,data(ii).uVel_Z,data(ii).vVel_Z,~] = CalcVorticity(MITfile,Tii);
    [data(ii).LON,data(ii).LAT] = ndgrid(data(ii).Lon,data(ii).Lat);
    dZLon = data(ii).ZLon(2,1)-data(ii).ZLon(1,1); dZLat = data(ii).ZLat(2,1)-data(ii).ZLat(1,1); 
    % barotropic streamfunction
    [data(ii).bsf,~] = CalcCavityStreamFunctions(runID,MITfile,Tii);

    load(Uafile);

    PVcutoff = -10e-7;
    data(ii).PV(data(ii).PV<PVcutoff)=PVcutoff;
    data(ii).PV(data(ii).PV>-1e-7)=-1e-7;
    contourf(ax_fig(ii),data(ii).ZLon/1e3,data(ii).ZLat/1e3,data(ii).PV,[PVcutoff:1e-7:-2e-7],'EdgeColor','none','LineWidth',0.5); 
%     CM1=[255,255,217;...
%         255,255,217;...
%         237,248,177;...
%         199,233,180;...
%         127,205,187;...
%         65,182,196;...
%         29,145,192;...
%         34,94,168;...
%         37,52,148;...
%         8,29,88]/255;
    %n = size(CM1,1)-1; 
    m = numel([PVcutoff:1e-7:-2e-7]);
%     CM = [interp1([0:n],CM1(:,1),[1:m]/m*n)',...
%         interp1([0:n],CM1(:,2),[1:m]/m*n)',...
%         interp1([0:n],CM1(:,3),[1:m]/m*n)'];
    CM = viridis_white(m);
    
    colormap(ax_fig(ii),CM);
    caxis(ax_fig(ii),[PVcutoff -1e-7]);

    Melt = double(ncread(MITfile,"SHIfwFlx"));
    Mask = Melt; Mask(Mask~=0)=1; Mask(Mask==0)=NaN;

    [cn,hn]=contour(ax_fig(ii),data(ii).LON/1e3,data(ii).LAT/1e3,data(ii).bsf.AS.fixedmask.*Mask/1e6,[-2:0.5:-0.5]*1e5/1e6,'--','Color',[0.2 0.2 0.2],'LineWidth',0.5);
    clabel(cn,hn,[-2:1:-0.5]*1e5/1e6,'labelspacing',200,'color','k','fontsize',7);
    [cz,hz]=contour(ax_fig(ii),data(ii).LON/1e3,data(ii).LAT/1e3,data(ii).bsf.AS.fixedmask.*Mask/1e6,[0 0],'-','Color',[0.2 0.2 0.2],'LineWidth',1);
    clabel(cz,hz,[0 0],'labelspacing',200,'color','k','fontsize',7);
    [cp,hp]=contour(ax_fig(ii),data(ii).LON/1e3,data(ii).LAT/1e3,data(ii).bsf.AS.fixedmask.*Mask/1e6,[0.5:0.5:8]*1e5/1e6,'-','Color',[0.2 0.2 0.2],'LineWidth',0.5);
    clabel(cp,hp,[0.5:1:8]*1e5/1e6,'labelspacing',200,'color','k','fontsize',7);

    nstr = strlength(MITfile);
    yyyymm = extractBetween(MITfile,nstr-22,nstr-17);
    time = datenum(yyyymm+"01","yyyymmdd");        

    [data(ii).xGL,data(ii).yGL,~]=PlotGroundingLines(CtrlVar,MUA,GF);
    data(ii).xB = MUA.Boundary.x;
    data(ii).yB = MUA.Boundary.y;

    [~,I] = min(abs(time-MITTime));
    data(ii).bsf.contours.xmid = section(I).xmid;
    data(ii).bsf.contours.ymid = section(I).ymid;

    plot(ax_fig(ii),data(ii).xB/1e3,data(ii).yB/1e3,'-k');
    plot(ax_fig(ii),data(ii).xGL/1e3,data(ii).yGL/1e3,'-k');
    plot(ax_fig(ii),[-1.605e6 -1.605e6 -1.58e6 -1.58e6 -1.605e6]/1e3,[-2.75e5 -2.5e5 -2.5e5 -2.75e5 -2.75e5]/1e3,':k','linewidth',3);
    plot(ax_fig(ii),data(ii).bsf.contours.xmid/1e3,data(ii).bsf.contours.ymid/1e3,'+m','markersize',4);

    axis(ax_fig(ii),'equal');
    xlim(ax_fig(ii),[xmin xmax]); ylim(ax_fig(ii),[ymin ymax]);
    grid(ax_fig(ii),'on'); box(ax_fig(ii),'on');

    title(ax_fig(ii),"year "+year(ii));

    if ii==1
        xticklabels(ax_fig(ii),'');
    elseif ii==2
        xticklabels(ax_fig(ii),'');
        yticklabels(ax_fig(ii),'');
    elseif ii==3
        xticklabels(ax_fig(ii),'');
        yticklabels(ax_fig(ii),'');
    elseif ii==4
        
    elseif ii==5
        yticklabels(ax_fig(ii),'');
    elseif ii==6
        yticklabels(ax_fig(ii),'');
    end

end

xlabel(tlo_fig,'psx [km]','fontsize',14); ylabel(tlo_fig,'psy [km]','fontsize',14);
cb = colorbar;
cb=colorbar("AxisLocation","out","Ticks",[PVcutoff:1e-7:-1e-7]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.02;
cb.FontSize=16;
cb.Label.String = "$\frac{f+\zeta}{H}\;\rm{[s}^{-1}\rm{m}^{-1}\rm{]}$";
cb.Label.Interpreter = "latex";
cb.Layout.Tile = 'east'; 


pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/FigS6_PIG_PV_and_StreamFunctions";   
print(H,fname,"-dpng","-r400");
