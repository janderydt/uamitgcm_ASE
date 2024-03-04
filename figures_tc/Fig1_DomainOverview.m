function Fig1_DomainOverview

%% Plotting routine for Fig 1 in https://doi.org/10.5194/egusphere-2023-1587
%% Requires the uamitgcm toolbox (https://github.com/janderydt/uamitgcm_tools)
%% and Ua Utilities (https://github.com/GHilmarG/UaSource/tree/beta/UaUtilities)

% Initialize UaMitgcm case directory and personal Antarctic EO data repository
froot = getenv("froot_uamitgcm");
froot_data = getenv("froot_data");

% Load UaMITgcm toolbox
addpath(getenv("froot_tools"));

%% Gather data
% Load MITgcm grids 
lonCMITbound = ncread(froot+"/cases/ASE_varmelt/output/221501/MITgcm/output.nc","XC")/1e3;
latCMITbound = ncread(froot+"/cases/ASE_varmelt/output/221501/MITgcm/output.nc","YC")/1e3;

% Load temperature data from Naughten et al. 2021 (https://doi.org/10.1029/2021GL094566)
lon=ncread('temp_max_below_400m_last10y_avg.nc','X');
lat=ncread('temp_max_below_400m_last10y_avg.nc','Y');
[Lonm,Latm] = ndgrid(lon,lat);
Xm=ncread('temp_max_below_400m_last10y_avg.nc','cartesian_x');
Ym=ncread('temp_max_below_400m_last10y_avg.nc','cartesian_y');
Theta_Naughten=ncread('temp_max_below_400m_last10y_avg.nc','THETA');
M=max(Theta_Naughten(:));
Theta_Naughten(Theta_Naughten==M)=NaN;
Theta_Naughten(Latm>-69.5)=NaN;

% Load IMBIE basin shapes
BasinShape = shaperead(froot_data+"/Antarctica_DrainageBasins_IMBIE/Basins_Antarctica_v02.shp");

% Load offshore bathymetry from BM Antarctica
addpath(froot_data+"/BedMachine_Antarctica/");
bathy_BMA = interpBedmachineAntarctica(Xm,Ym,'bed',"2020-07-15");
bathy_BMA = reshape(bathy_BMA,size(Xm));
Ibma = find(Latm>-69.5); 
bathy_BMA(Ibma)=nan;

% Obtain Antarctic coastline from Bedmap2
addpath(genpath(froot_data+"/Bedmap2/"));
[x,y,~] = bedmap2_data('surface','xy');
mask = bedmap2_data('icemask');
[C,h] = contour(x,y,mask,[127 127]);
ii=1; kk=1;
while ii<size(C,2)
   step = C(2,ii);
   Cont(kk).step = step;
   Cont(kk).x = C(1,ii+1:ii+step);
   Cont(kk).y = C(2,ii+1:ii+step);
   ii = ii + step + 1;
   kk = kk+1;   
end
[Maxstepsize,I] = max([Cont(:).step]);
ContAnt_x = Cont(I).x(:);
ContAnt_y = Cont(I).y(:);

% Remove bathymetry values for ice-covered regions
I = find(inpoly([Xm(:),Ym(:)],[ContAnt_x(:) ContAnt_y(:)]));
bathy_BMA(I)=nan;

% Load Ua data
load("./UaDefaultRun_199701_0031.mat");
MUA.coordinates = MUA.coordinates/1e3;

CtrlVar.PlotNodes=0;
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
CtrlVar.MeshColor = [0.6 0.6 0.6];
CtrlVar.FEmeshPlotTitle="";
CtrlVar.PlotsXaxisLabel = "psx [km]";
CtrlVar.PlotsYaxisLabel = "psy [km]";
CtrlVar.PlotXYscale = 1;

%% Start plotting
% Set bounds for the figure 
xmin = -2.25e6; xmax = -0.99e6; 
ymin = -1.1e6; ymax = 0.2e6;

% New figure
H=fig('units','inches','width',60*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');

subplot("position",[0.03 0.05 0.85 0.92]); hold on;

% Ocean temperature 
alphaVal = 0.5;
pcolor(Xm/1e3,Ym/1e3,Theta_Naughten); shading flat; alpha(alphaVal);
caxis(gca,[-0.5 1.5]);

% Bathymetric contours
[C,h]=contour(Xm/1e3,Ym/1e3,bathy_BMA,[-4500:500:-1000 -400],"-","color",[0.5 0.5 0.5]);
clabel(C,h,[-4000:1000:-1000 -400],"labelspacing",400,"color",[0.5 0.5 0.5],"fontsize",8);

% LatLon Grid (script adapted from Ua source code https://github.com/GHilmarG/UaSource/)
PlotLatLonGrid(1e3, 2.5, 15, 1.2, [0.3 0.3 0.3], 0.5, 11);

% Plot Ua Mesh (requires Ua source code https://github.com/GHilmarG/UaSource/)
PlotMuaMesh(CtrlVar,MUA);

% Plot Ua and MITgcm boundaries
g(1) = plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"k","linewidth",2);
g(3) = plot([-1.6e3  min(lonCMITbound(:)) min(lonCMITbound(:))],...
    [min(latCMITbound(:)) min(latCMITbound(:)) -3.4e2],...
    "color",[0.31 0.28 0.7],"linestyle",":","linewidth",4);%[246,190,0]/255
g(2) = plot([min(lonCMITbound(:)) min(lonCMITbound(:)) max(lonCMITbound(:)) max(lonCMITbound(:)) -1.6e3],...
    [-3.4e2 max(latCMITbound(:)) max(latCMITbound(:)) min(latCMITbound(:)) min(latCMITbound(:))],...
    "color",[0.31 0.28 0.7],"linestyle","-","linewidth",4);%[246,190,0]/255

% Plot Ua grounding lines (requires Ua source code https://github.com/GHilmarG/UaSource/)
PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"k","linewidth",2);

% Plot IMBIE basin shapes
for ii=1:length(BasinShape)
     X = double(BasinShape(ii).X)/1e3;
     Y = double(BasinShape(ii).Y)/1e3;
    plot(X,Y,"color",[0.4 0.4 0.4]); 
end

I = find(inpoly([ContAnt_x(:)/1e3 ContAnt_y(:)/1e3],[[min(lonCMITbound(:)) min(lonCMITbound(:)) max(lonCMITbound(:)) max(lonCMITbound(:))]',...
    [min(latCMITbound(:)) max(latCMITbound(:)) max(latCMITbound(:)) min(latCMITbound(:))]']));
ContAnt_x(I)=NaN; ContAnt_y(I)=NaN;
plot(ContAnt_x/1e3,ContAnt_y/1e3,"-k","linewidth",1);

patch([xmin xmin xmax xmax]/1e3,...
    [ymin ymin+200e3 ymin+200e3 ymin]/1e3,"w","edgecolor","w");

% Plot length scale - requires UaMITgcm toolbox
PlotLengthScale(xmin/1e3+900,ymin/1e3+230,[0 100 200 300],10);

% Deal with colorbars, titles and axes
dCM = 64;
CM2 = othercolor("BuOrR_14",2*dCM); CM2 = CM2(dCM+5:end,:);
dCM = 16;
CM1 = [linspace(1,CM2(1,1),dCM+1)' linspace(1,CM2(1,2),dCM+1)' linspace(1,CM2(1,3),dCM+1)'];
dCM = 16;
CM3 = [linspace(CM2(end,1),0.7,dCM+1)' linspace(CM2(end,2),0.7,dCM+1)' linspace(CM2(end,3),0.7,dCM+1)'];
CM = [CM1;CM2;CM3(2:end,:)];
CM2 = CM;
CM2 = min(1.3*CM2,zeros(size(CM2))+1);
colormap(gca,CM2);
cb=colorbar("Position",[0.125 0.254 0.012 0.25],"Location","west","AxisLocation","out","Ticks",[-0.5:0.5:1.5]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.04;
cb.FontSize=14;
cb.Label.String = "$\mbox{Maximum temperature below -400m [}^{\circ}\mbox{C]}$";
cb.Label.Interpreter = "latex";

title(""); xlabel(gca,""); ylabel(gca,"");
set(gca,"xtick",[],"ytick",[],"xticklabels",{},"yticklabels",{});

xlim([xmin xmax]/1e3); ylim([ymin ymax]/1e3);

box off; grid off;

legend(g(:),{"$\mbox{{\'U}a ice boundary}$","$\mbox{MITgcm ocean boundary}$","$\mbox{MITgcm open boundary}$"},"location","northeast","interpreter","latex","fontsize",14);

% Save
pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/Figure1";
print(H,fname,"-dpng","-r400");

end