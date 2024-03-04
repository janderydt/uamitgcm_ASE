function FigA2_UaValidation

runID="ASE_varmelt";
basin = "AS";
xmin = -1750; xmax = -1350;
ymin = -750; ymax = 0;

addpath(getenv("froot_tools"));
froot_UaMITgcm = getenv("froot_uamitgcm");
froot_data = getenv("froot_data");
frootm = froot_UaMITgcm+"/cases/"+runID;

timestep1 = 1;%*12;
timestep2 = 68*12;%205*12;
timestep3 = 118*12;

if timestep1>=timestep2
   error("check timesteps"); 
end

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = {subd(isub).name}';
nameFolds(ismember(nameFolds,[".",".."])) = [];

restartfile = dir(frootm+"/ua_custom/*-RestartFile.mat");
files{1}.folder = restartfile(1).folder;
files{1}.name = restartfile(1).name;

nfolders = length(nameFolds); kk=2;
for jj=1:timestep2
    %disp([frootm,"/output/",nameFolds{jj}]);
    output(jj).filelist=dir(frootm+"/output/"+nameFolds{jj}+"/Ua/UaDefaultRun_"+char(nameFolds{jj})+"*.mat");
    for ii=1:length(output(jj).filelist)
        files{kk}.folder = output(jj).filelist(ii).folder;
        files{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
end

load(files{timestep1}.folder+"/"+files{timestep1}.name);
year0 = files{timestep1+1}.name(end-14:end-11);
month0 = files{timestep1+1}.name(end-10:end-9);
if timestep1==1
    CtrlVar = CtrlVarInRestartFile;
    UserVar = UserVarInRestartFile;
    s = F.s; ub = F.ub; vb = F.vb; b = F.b; B=F.B;
    load(files{timestep1+1}.folder+"/"+files{timestep1+1}.name,"time");
end
time0 = datestr(datenum("01/"+month0+"/"+year0,"dd/mm/yyyy"),"dd/mm/yyyy")
    
MUA0 = MUA; GF0=GF;
CtrlVar.PlotXYscale = 1e3;
Fs0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),s);
speed0 = sqrt(ub.^2+vb.^2);
Fspeed0 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),speed0);

if timestep2 > length(files)
    timestep2 = length(files);
end

year1 = files{timestep2}.name(end-14:end-11);
month1 = files{timestep2}.name(end-10:end-9);
days1 = str2double(files{timestep2}.name(end-7:end-4));
datenum("01/"+month1+"/"+year1);
time1 = datestr(datenum("01/"+month1+"/"+year1,"dd/mm/yyyy")+days1,"dd/mm/yyyy")
load(files{timestep2}.folder+"/"+files{timestep2}.name,"h","s","ub","vb","MUA","GF","ab","b","B");
%deltah = h-Fh0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltas = s-Fs0(MUA.coordinates(:,1),MUA.coordinates(:,2));
deltaspeed = sqrt(ub.^2+vb.^2)-Fspeed0(MUA.coordinates(:,1),MUA.coordinates(:,2));

ub(h<=1.5)=NaN; vb(h<=1.5)=NaN;

xmin2 = min(MUA.coordinates(:,1)); xmax2 = max(MUA.coordinates(:,1));
ymin2 = min(MUA.coordinates(:,2)); ymax2 = max(MUA.coordinates(:,2));

CM = othercolor("RdBu11");
CM1 = othercolor("PuRd9",8); CM1 = flipdim(CM1,1);
CM2 = othercolor("Blues9",8); %
CM = [CM1 ; CM2]; CM = flipdim(CM,1);
%CM=flipdim(CM,1);


%% plotting 
Hfig = fig("units","inches","height",47*12/72.27,"width",90*12/72.27,"fontsize",12,"font","Helvetica");

tlo_fig = tiledlayout(1,4,"TileSpacing","compact");
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

%% Panel 1

PlotNodalBasedQuantities_JDR(ax_fig(1),MUA.connectivity,MUA.coordinates/1e3,-deltas); hold on;

colormap(ax_fig(1),CM); caxis(ax_fig(1),[-50 50]); 
            
plot(ax_fig(1),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");
CtrlVar.PlotGLs=0;
[xGL0,yGL0,~]=PlotGroundingLines(CtrlVar,MUA0,GF0);
g1=plot(ax_fig(1),xGL0/1e3,yGL0/1e3,"k","linewidth",1.5);
[xGL,yGL,~]=PlotGroundingLines(CtrlVar,MUA,GF);
g2=plot(ax_fig(1),xGL/1e3,yGL/1e3,"c","linewidth",1);
%PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"m","linewidth",1.5);

legend(ax_fig(1),[g1 g2],["$\mbox{{\'U}aMITgcm 1997}$","$\mbox{{\'U}aMITgcm 2015}$"],...
    "location","northeast","fontsize",12,"Interpreter","latex");

grid(ax_fig(1),"on"); box(ax_fig(1),"on");
axis(ax_fig(1),"equal");
xlim(ax_fig(1),[xmin xmax]); ylim(ax_fig(1),[ymin ymax]);

%xticklabels(ax_fig(1),"");
title(ax_fig(1),{"$\mbox{Simulated (\'UaMITgcm)}$";"1997-2014"},"interpreter","latex","fontsize",14);

%% Panel 2
%% CPOM
addpath(froot_data+"/CPOM_dhdt");
[dh_CPOM] = Interpolate_CPOMdhdt(MUA,["1997_2001", "2002_2006", "2007_2011", "2012_2016"]);

dh_CPOM(GF.node<1)=NaN;
PlotNodalBasedQuantities_JDR(ax_fig(2),MUA.connectivity,MUA.coordinates/1e3,-dh_CPOM); hold on;

%% ITS_LIVE

cb1=colorbar(ax_fig(1),"Ticks",[-50:25:50],"Position",[0.13 0.08 0.3 0.01],"Location","south","AxisLocation","out");
cb1.XColor="k";
cb1.YColor="k";
cb1.TickLength=0.02;
cb1.FontSize=12;
cb1.Label.String = "Surface elevation change [m]";

colormap(ax_fig(2),CM); caxis(ax_fig(2),[-50 50]);

plot(ax_fig(2),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");
%PlotGroundingLines(CtrlVar,MUA0,GF0,[],[],[],"k","linewidth",1.5);
%PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"b","linewidth",1.5);

S = shaperead(froot_data+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
for ii=1:length(S)
   date=S(ii).DATE1;
   year = date(1:4);
   [x,y] = ll2psxy(S(ii).Y,S(ii).X,-71,0);
   if ismember(str2double(year),[1992:1996])
        g3=plot(ax_fig(2),x/1e3,y/1e3,"k","linewidth",1.5);
   end
   if str2double(year)==2011
        g4=plot(ax_fig(2),x/1e3,y/1e3,"c","linewidth",1.5); 
   end
end

legend(ax_fig(2),[g3 g4],["$\mbox{InSAR 1992-1996}$","$\mbox{InSAR 2012}$"],...
    "location","northeast","fontsize",12,"interpreter","latex");

grid(ax_fig(2),"on"); box(ax_fig(2),"on"); axis equal;
axis(ax_fig(2),"equal");
xlim(ax_fig(2),[xmin xmax]); ylim(ax_fig(2),[ymin ymax]);

%xticklabels(ax_fig(2),""); 
yticklabels(ax_fig(2),"");
title(ax_fig(2),{"$\mbox{Measured (CPOM)}$";"1997-2016"},"interpreter","latex","fontsize",14);
  
%% Panel 3

PlotNodalBasedQuantities_JDR(ax_fig(3),MUA.connectivity,MUA.coordinates/1e3,deltaspeed);

colormap(ax_fig(3),CM); caxis(ax_fig(3),[-500 500]);
 

%              N=10;
%              quiver(Xm(1:N:end,1:N:end)/1e3,Ym(1:N:end,1:N:end)/1e3,du(1:N:end,1:N:end)/1e2,dv(1:N:end,1:N:end)/1e2,"k","autoscale","off","linewidth",1);
%              quiver(xmin+10,ymax-15,1000/1e2,0,"k","autoscale","off","linewidth",2,"MaxHeadSize",1);
%              text(xmin+10,ymax-30,"1000 m/yr");
        
   
plot(ax_fig(3),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");

plot(ax_fig(3),xGL0/1e3,yGL0/1e3,"k","linewidth",1.5);
plot(ax_fig(3),xGL/1e3,yGL/1e3,"c","linewidth",1);

grid(ax_fig(3),"on"); box(ax_fig(3),"on");
axis(ax_fig(3),"equal");
xlim(ax_fig(3),[xmin xmax]); ylim(ax_fig(3),[ymin ymax]);
yticklabels(ax_fig(3),"");

title(ax_fig(3),{"$\mbox{Simulated (\'UaMITgcm)}$";"1997-2014"},"interpreter","latex","fontsize",14);

%% Panel 4

addpath(froot_data+"/Measures/");
addpath(froot_data+"/PIG_Thwaites_DotsonCrosson/VelocityMaps/ALOS, ERS-1, RADARSAT-1, RADARSAT-2, TDX");
data = ReadNetCDF_ASE;
X1=data(1).psx; Y1=data(1).psy; VX1=data(1).vx; VY1=data(1).vy;
[X1,Y1]=ndgrid(X1,Y1);
Fv1996 = griddedInterpolant(X1,Y1,sqrt(VX1.^2+VY1.^2));
[X2,Y2,VX2,VY2,~] = ReadMeasures("2014_2015",[xmin xmax ymin ymax]*1e3);
Fv2015 =  griddedInterpolant(X2,Y2,sqrt(VX2.^2+VY2.^2));
dspeed_Meas = Fv2015(MUA.coordinates(:,1),MUA.coordinates(:,2))-Fv1996(MUA.coordinates(:,1),MUA.coordinates(:,2));

%dspeed_Meas(GF.node<1)=NaN;
PlotNodalBasedQuantities_JDR(ax_fig(4),MUA.connectivity,MUA.coordinates/1e3,dspeed_Meas); hold on;

colormap(ax_fig(4),CM); caxis(ax_fig(4),[-500 500]);

cb2=colorbar(ax_fig(3),"Ticks",[-500:250:500],"Position",[0.55 0.08 0.3 0.01],"Location","south","AxisLocation","out");
cb2.XColor="k";
cb2.YColor="k";
cb2.TickLength=0.02;
cb2.FontSize=12;
cb2.Label.String = "Surface velocity change [m yr^{-1}]";

plot(ax_fig(4),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k");
S = shaperead(froot_data+"GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
for ii=1:length(S)
   date=S(ii).DATE1;
   year = date(1:4);
   [x,y] = ll2psxy(S(ii).Y,S(ii).X,-71,0);
   if ismember(str2double(year),[1992:1996])
        plot(ax_fig(4),x/1e3,y/1e3,"k","linewidth",1.5);
   end
   if str2double(year)==2011
        plot(ax_fig(4),x/1e3,y/1e3,"c","linewidth",1.5); 
   end
end

grid(ax_fig(4),"on"); box(ax_fig(4),"on");
axis(ax_fig(4),"equal");
xlim(ax_fig(4),[xmin xmax]); ylim(ax_fig(4),[ymin ymax]);
yticklabels(ax_fig(4),"");

title(ax_fig(4),{"$\mbox{Measured (MeaSUREs)}$";"1995-2015"},"interpreter","latex","fontsize",14);

xlabel(tlo_fig,"psx [km]","Fontsize",12);
ylabel(tlo_fig,"psy [km]","Fontsize",12);

%Pos = cb1.Position; cb1.Position=[Pos(1)+0.03 Pos(2) Pos(3)-0.008 Pos(4)];
%Pos = cb2.Position; cb2.Position=[Pos(1)+0.03 Pos(2) Pos(3)-0.008 Pos(4)];

%% Save
pos = get(Hfig,"Position");
set(Hfig,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);

folder = "../Figures/";
if ~exist(folder)
    mkdir(folder);
end
fname = folder+"/FigA2_Change_in_SurfaceEl_Speed"+"_"+basin;
print(Hfig,fname,"-dpng","-r400");