function Fig4_FigS1_ChangeInWCT

%% Plotting routine for Figs 4 and S1 in https://doi.org/10.5194/egusphere-2023-1587
%% Requires the uamitgcm toolbox (https://github.com/janderydt/uamitgcm_tools)
%% and Ua Utilities (https://github.com/GHilmarG/UaSource/tree/beta/UaUtilities)


runID="ASE_himelt";
basins = ["PIG","TW","DC"];

panellabels = {'a','b','c'};
timestep1 = 2; % after 5-year relaxation
timestep2 = 155*12+1;

addpath(getenv("froot_tools"));
froot_data = getenv("froot_data");
frootm = getenv("froot_uamitgcm")+"cases/"+runID;
%frootm = ["/media/janryd69/mainJDeRydt/UaMITgcm_v2/cases/",runID];

if timestep1>=timestep2
   error("check timesteps"); 
end

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = {subd(isub).name}';
nameFolds(ismember(nameFolds,[".",".."])) = [];
nfolders = numel(nameFolds); kk=1; jj=1;
while kk<=timestep2 && jj<nfolders
    disp(frootm+"/output/"+nameFolds{jj});
    output(jj).filelist=dir(frootm+"/output/"+nameFolds{jj}+"/Ua/UaDefaultRun_"+char(nameFolds{jj})+"*.mat");
    for ii=1:length(output(jj).filelist)
        files{kk}.folder = output(jj).filelist(ii).folder;
        files{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
    jj=jj+1;
end

restartfile = dir(frootm+"/ua_custom/*-RestartFile.mat");
files{1}.folder = restartfile(1).folder;
files{1}.name = restartfile(1).name;
%[files{timestep1}.folder,"/",files{timestep1}.name]
load(files{timestep1}.folder+"/"+files{timestep1}.name);
year0 = files{timestep1+1}.name(end-14:end-11);
month0 = files{timestep1+1}.name(end-10:end-9);
if timestep1==1
    CtrlVar = CtrlVarInRestartFile;
    UserVar = UserVarInRestartFile;
    rho = F.rho; s = F.s; h = F.h; ub = F.ub; vb = F.vb; ab=F.ab;
    b = F.b; B=F.B;
    load(files{timestep1+1}.folder+"/"+files{timestep1+1}.name+"time");
end
time0 = datestr(datenum("01/"+month0+"/"+year0,"dd/mm/yyyy"),"dd/mm/yyyy");
    
MUA0 = MUA; GF0=GF;
CtrlVar.PlotXYscale = 1e3;
Fh0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),h);
Fs0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),s);
Fwct0 =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),b-B);
speed0 = sqrt(ub.^2+vb.^2);
Fspeed0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed0);
Fub0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ub);
Fvb0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),vb);
Fab0 =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ab);

%if timestep2 > length(files)
%    timestep2 = length(files);
%end
files{timestep2}.folder,"/",files{timestep2}.name
year1 = files{timestep2}.name(end-14:end-11);
month1 = files{timestep2}.name(end-10:end-9);
days1 = str2double(files{timestep2}.name(end-7:end-4));
datenum("01/"+month1+"/"+year1);
time1 = datestr(datenum("01/"+month1+"/"+year1,"dd/mm/yyyy")+days1,"dd/mm/yyyy");
load(files{timestep2}.folder+"/"+files{timestep2}.name,"h","s","ub","vb","MUA","GF","ab","b","B");
%deltah = h-Fh0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltah = s-Fs0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltaspeed = sqrt(ub.^2+vb.^2)-Fspeed0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltaab = ab-Fab0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltawct = b-B-Fwct0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltawct(GF.node==1)=NaN;

Fspeed =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),sqrt(ub.^2+vb.^2));
Bcont = B; Bcont(GF.node<1)=NaN;
FBcont =  scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),Bcont);
wctcont = b-B; wctcont(GF.node==1)=NaN;
Fwctcont = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),wctcont);
Fdeltaspeed = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),deltaspeed);
ub(h<=1.5)=NaN; vb(h<=1.5)=NaN;
Fub = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ub);
Fvb = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),vb);
Fdeltah = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),deltah);

xmin2 = min(MUA.coordinates(:,1)); xmax2 = max(MUA.coordinates(:,1));
ymin2 = min(MUA.coordinates(:,2)); ymax2 = max(MUA.coordinates(:,2));

[Xm,Ym] = ndgrid([xmin2:500:xmax2],[ymin2:500:ymax2]);
deltaspeedm = Fdeltaspeed(Xm,Ym);
deltahm = Fdeltah(Xm,Ym); 
speedm = Fspeed(Xm,Ym);
Bm = FBcont(Xm,Ym);
wctm = Fwctcont(Xm,Ym);
J = find(isnan(MUA.Boundary.x));
if isempty(J)
    J=length(MUA.Boundary.x)+1;
end
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x(1:J(1)-1) MUA.Boundary.y(1:J(1)-1)]));
K = find(abs(deltaspeedm) < 20);
deltaspeedm(I) = NaN;
deltahm(I)=NaN;
du = Fub(Xm,Ym)-Fub0(Xm,Ym); du([I;K])= NaN;
dv = Fvb(Xm,Ym)-Fvb0(Xm,Ym); dv([I;K])= NaN;
speedm(I) = NaN;
Bm(I) = NaN;
wctm(I) = NaN;


S_insar = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
S_insar_mililloTW = shaperead(froot_data+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
S_insar_mililloDC = shaperead(froot_data+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");


%% plotting 
Hfig = fig('units','inches','height',42.84*12/72.27,'width',107.1*12/72.27,'fontsize',12,'font','Helvetica');

tlo_fig = tiledlayout(1,3,"TileSpacing","tight");
for i = 1:3
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

for bb=1:numel(basins)

    switch basins{bb}
        case "PIG"
            %% PIG
            xmin = -1700; xmax = -1530; 
            ymin = -370; ymax = -210;
            basinname="Pine Island";
        case "TW"
            %% Thwaites
            xmin = -1620; xmax = -1450;
            ymin = -530; ymax = -370;
            basinname = "Thwaites";
        case "DC"
            %% Dotson Crosson
            xmin = -1620; xmax = -1440;
            ymin = -700; ymax = -530;
            basinname = "Crosson - Dotson";
        case "AS"
            %% Amundsen
            xmin = -1700; xmax = -1440;
            ymin = -700; ymax = -210;
        otherwise
            error("case unknown.")
    end

            
            
    PlotNodalBasedQuantities_JDR(ax_fig(bb),MUA.connectivity,MUA.coordinates/1e3,deltawct); hold on;

    Bcont(Bcont<-1600)=-1600; Bcont(Bcont>-200)=-200;
    scale = (-Bcont+max(Bcont,[],"omitnan"))/(max(Bcont,[],"omitnan")-min(Bcont,[],"omitnan"));
    Mask = GF.node;
    alphaMask = GF.node; I = find(alphaMask<1);
    alphaMask(I)=0;

    patch(ax_fig(bb),"faces",MUA.connectivity,"vertices",MUA.coordinates/1e3,...
    "FaceVertexCData",Mask,"CDataMapping","scaled","EdgeColor","none","FaceColor",[0.5 0.5 0.8],...
    "FaceAlpha","flat","FaceVertexAlphaData",alphaMask.*scale);

    [M2,c2]=contour(ax_fig(bb),Xm/1e3,Ym/1e3,Bm,[-1600:100:-600]); c2.Color=[0.5 0.5 1];
    clabel(M2,c2,[-1600:200:-200],"color",[0.3 0.3 1],"labelspacing",500); 

    
    for ii=1:length(S_insar)
        date=S_insar(ii).DATE1;
        year_insar = date(1:4);
        [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
        if str2double(year_insar)>=2011
            ginsar=plot(ax_fig(bb),x_insar/1e3,y_insar/1e3,"color",[255, 165, 0]/255,"linewidth",2.5); 
        end
    end

    for ii=1:length(S_insar_mililloTW)
        plot(ax_fig(bb),S_insar_mililloTW(ii).X/1e3,S_insar_mililloTW(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);  
    end

    for ii=1:length(S_insar_mililloDC)
        plot(ax_fig(bb),S_insar_mililloDC(ii).X/1e3,S_insar_mililloDC(ii).Y/1e3,"color",[255, 165, 0]/255,"linewidth",2.5);                
    end
    
    plot(ax_fig(bb),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");
    
    CtrlVar.PlotGLs=0; CtrlVar.PlotXYscale = 1e3;
    [xGL,yGL,~]=PlotGroundingLines(CtrlVar,MUA,GF);
    gl_200 = plot(ax_fig(bb),xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,"m","linewidth",1.5) ;
    [xGL0,yGL0,~]=PlotGroundingLines(CtrlVar,MUA0,GF0);
    gl_0 = plot(ax_fig(bb),xGL0/CtrlVar.PlotXYscale,yGL0/CtrlVar.PlotXYscale,"k","linewidth",1.5) ;

    if bb==1
        plot(ax_fig(bb),1e3*[-1.5859  -1.619],1e3*[-0.240 -0.3686],'--k','LineWidth',2);
    elseif bb==3
        plot(ax_fig(bb),1e3*[-1.475 -1.495  -1.5306],1e3*[-0.662 -0.67 -0.6529],'--k','LineWidth',2);
    end

    if bb==2
        legend(ax_fig(bb),[ginsar,gl_0,gl_200],{'$\mbox{InSAR 2011-2017}$','$\mbox{{\''U}a-MITgcm year 0}$','$\mbox{{\''U}a-MITgcm year 200}$'},"location","northeast","fontsize",12,"interpreter","latex");
    end

    %text(ax_fig(bb),xmin+5e3,ymax-5e3,panellabels{bb},'HorizontalAlignment','left','VerticalAlignment','top','fontsize',14);

    axis(ax_fig(bb),"equal"); box(ax_fig(bb),"on");
    xlim(ax_fig(bb),[xmin xmax]); ylim(ax_fig(bb),[ymin ymax]);



    %set([ax1],"Position",[0.1 0.1 0.85 0.85]);

    %% Give each one its own colormap
    CM = othercolor("RdYlBu11",23); CM = CM(12:23,:); CM(2,:)=[1 1 1];
    colormap(ax_fig(bb),CM); caxis(ax_fig(bb),[-100 750]);

    if contains(basins{bb},"PIG")
        cb1=colorbar(ax_fig(bb),"Position",[0.1 0.79 0.10 0.015],"Location","northoutside","AxisLocation","in","Ticks",[-100 0:250:750]);
        cb1.XColor="k";
        cb1.YColor="k";
        cb1.TickLength=0.04;
        cb1.FontSize=12;
        cb1.Label.String = "Water thickness change [m]";
    end

    %CM2 = othercolor("Bu_10",6); 
    %colormap(ax2,CM2); caxis(ax2,[-1600 -600]);
    
%            title("UaMITgcm");
    if bb==1
            ylabel(ax_fig(bb),"psy [km]","Fontsize",12);
    elseif bb==2
            xlabel(ax_fig(bb),"psx [km]","Fontsize",12);
    end

    title(ax_fig(bb),basinname);

    end
            
end
    
pos = get(Hfig,"Position");
set(Hfig,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);

folder = "../Figures/"+runID;
if ~exist(folder)
    mkdir(folder);
end
if strcmp(runID,"ASE_himelt")
    fname = folder+"/Fig4_ChangeInWCT_ASE_himelt";
elseif strcmp(runID,"ASE_avmelt")
    fname = folder+"/FigS2_ChangeInWCT_ASE_avmelt";
end
print(Hfig,fname,"-dpng","-r400");