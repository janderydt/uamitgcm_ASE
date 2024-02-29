function Fig1_DomainOverview

froot = getenv("froot_uamitgcm");
froot_data = getenv("froot_data");

%lonCMIT = ncread(froot+"/cases/PTDC_001/output/200001/MITgcm/output.nc","XC")/1e3;
%latCMIT = ncread(froot+"/cases/PTDC_001/output/200001/MITgcm/output.nc","YC")/1e3;
lonCMIT = ncread(froot+"/cases/PTDC_002/output/275001/MITgcm/output.nc","XC")/1e3;
latCMIT = ncread(froot+"/cases/PTDC_002/output/275001/MITgcm/output.nc","YC")/1e3;
[lonCMIT,latCMIT] = ndgrid(lonCMIT,latCMIT);
hFacC = ncread(froot+"/cases/PTDC_001/output/200001/MITgcm/output.nc","hFacC");
%xmin = -2.5e6; xmax = -0.5e6; 
%ymin = -1.5e6; ymax = 1.7e6;
xmin = -2.25e6; xmax = -0.99e6; 
ymin = -1.1e6; ymax = 0.2e6;

addpath(getenv("froot_tools"));
Basins = DefineBasins;

%% B-SOSE data
% file_sose = froot_data+"/B-SOSE/bsose_i105_2008to2012_monthly_Theta.nc";
% lonC_sose = ncread(file_sose,"XC");
% latC_sose = ncread(file_sose,"YC");
% [LONC_sose,LATC_sose] = ndgrid(lonC_sose,latC_sose);
% I69deg = find(LATC_sose>-69.5); 
% [XC_sose,YC_sose] = ll2psxy(LATC_sose,LONC_sose,-71,0);
% [Xm,Ym] = ndgrid(linspace(xmin,xmax,500),linspace(ymin,ymax,500));
% [Latm,Lonm] = psxy2ll(Xm,Ym,-71,0);
% % 
% depth_sose = ncread(file_sose,"Z");
% I = min(find(depth_sose<=-400));
% 
% Theta_sose = ncread(file_sose,"THETA");
% % find maximum average temperature below -400m
% Theta_sose_av = mean(Theta_sose(:,:,I:end,:),4);
% Theta_sose_av_max = max(Theta_sose_av,[],3);
% Theta_sose_av_max(I69deg) = NaN;
% FTheta = scatteredInterpolant(XC_sose(:),YC_sose(:),Theta_sose_av_max(:));
% Theta_sose = FTheta(Xm,Ym);

lon=ncread('temp_max_below_400m_last10y_avg.nc','X');
lat=ncread('temp_max_below_400m_last10y_avg.nc','Y');
[Lonm,Latm] = ndgrid(lon,lat);
Xm=ncread('temp_max_below_400m_last10y_avg.nc','cartesian_x');
Ym=ncread('temp_max_below_400m_last10y_avg.nc','cartesian_y');
Theta_Naughten=ncread('temp_max_below_400m_last10y_avg.nc','THETA');
M=max(Theta_Naughten(:));
Theta_Naughten(Theta_Naughten==M)=NaN;
Theta_Naughten(Latm>-69.5)=NaN;

figure; pcolor(Xm,Ym,Theta_Naughten); shading flat;

%% offshore bathymetry
addpath(froot_data+"/BedMachine_Antarctica/");
bathy_BMA = interpBedmachineAntarctica(Xm,Ym,'bed',"2020-07-15");
bathy_BMA = reshape(bathy_BMA,size(Xm));
Ibma = find(Latm>-69.5);
%bathy_BMA = flipdim(bathy_BMA,2); 
bathy_BMA(Ibma)=NaN;
figure; pcolor(Xm,Ym,bathy_BMA); shading flat;

%% Antarctic outline
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

I = find(inpoly([Xm(:),Ym(:)],[ContAnt_x(:) ContAnt_y(:)]));
bathy_BMA(I)=NaN;
Theta_sose(I)=NaN;

%% Toshi model domain
% Lat_T = [-75.5*ones(1,100) linspace(-75.5,-62,100) -62*ones(1,100) linspace(-62,-75.5,100)];
% Lon_T = [-linspace(140,80,100) -80*ones(1,100) -linspace(80,140,100) -140*ones(1,100)];
% [X_T,Y_T] = ll2psxy(Lat_T,Lon_T,-71,0);
% I = find(inpoly([X_T(:),Y_T(:)],[ContAnt_x(:) ContAnt_y(:)]));
% X_T(I)=NaN; Y_T(I)=NaN;

hFacC_vertint = sum(hFacC,3);
mask = 0*hFacC_vertint;
mask(hFacC_vertint==0)=0; %grounded
mask(hFacC_vertint>0)=1; % ice shelf
mask(hFacC(:,:,1)==1)=2; % open ocean

load("./UaDefaultRun_199701_0031.mat");
MUA.coordinates = MUA.coordinates/1e3;

CtrlVar.PlotNodes=0;
CtrlVar.WhenPlottingMesh_PlotMeshBoundaryCoordinatesToo=0;
CtrlVar.MeshColor = [0.6 0.6 0.6];
CtrlVar.FEmeshPlotTitle="";
CtrlVar.PlotsXaxisLabel = "psx [km]";
CtrlVar.PlotsYaxisLabel = "psy [km]";
CtrlVar.PlotXYscale = 1;

lonCMIT_orig = lonCMIT; latCMIT_orig = latCMIT;
lonCMIT(mask==0)=NaN;
latCMIT(mask==0)=NaN;
   
H=fig('units','inches','width',60*12/72.27,'height',60*12/72.27,'fontsize',14,'font','Helvetica');
%set(H,"visible","off");

subplot("position",[0.03 0.05 0.85 0.92]); hold on;

%lima(-1400e3,-350e3,4000,"xy");
alphaVal = 0.5;
%pcolor(Xm/1e3,Ym/1e3,Theta_sose); shading flat; alpha(alphaVal);
pcolor(Xm/1e3,Ym/1e3,Theta_Naughten); shading flat; alpha(alphaVal);
%caxis([min(Theta_sose(:)) max(Theta_sose(:))]);
CM = othercolor("RdBu6",64); CM = flipdim(CM,1); CM = CM(15:end,:);

%CM1 = othercolor("BuOrR_14",2*dCM); CM1 = CM1(1:dCM-5,:);

dCM = 64;
CM2 = othercolor("BuOrR_14",2*dCM); CM2 = CM2(dCM+5:end,:);
dCM = 16;
CM1 = [linspace(1,CM2(1,1),dCM+1)' linspace(1,CM2(1,2),dCM+1)' linspace(1,CM2(1,3),dCM+1)'];
dCM = 16;
%CM3 = [linspace(CM2(end,1),152/255,dCM+1)' linspace(CM2(end,2),68/255,dCM+1)' linspace(CM2(end,3),158/255,dCM+1)'];
CM3 = [linspace(CM2(end,1),0.7,dCM+1)' linspace(CM2(end,2),0.7,dCM+1)' linspace(CM2(end,3),0.7,dCM+1)'];
CM = [CM1;CM2;CM3(2:end,:)];

caxis([-0.5 1.5]);

colormap(gca,CM);

[C,h]=contour(Xm/1e3,Ym/1e3,bathy_BMA,[-4500:500:-1000 -400],"-","color",[0.5 0.5 0.5]);
clabel(C,h,[-4000:1000:-1000 -400],"labelspacing",400,"color",[0.5 0.5 0.5],"fontsize",8);

%plot(X_T/1e3,Y_T/1e3,"--k","linewidth",2);

PlotLatLonGrid(1e3, 2.5, 15, 1.2, [0.3 0.3 0.3], 0.5, 11);

%plot(lonCMIT,latCMIT,"color",[246,190,0]/255);
%plot(lonCMIT',latCMIT',"color",[246,190,0]/255);
PlotMuaMesh(CtrlVar,MUA);

g(1) = plot(MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"k","linewidth",2);

lonCMITbound = ncread(froot+"/cases/PTDC_001/output/221501/MITgcm/output.nc","XC")/1e3;
latCMITbound = ncread(froot+"/cases/PTDC_001/output/221501/MITgcm/output.nc","YC")/1e3;
%lonCMITbound = ncread(froot+"/cases/PTDC_002/output/275001/MITgcm/output.nc","XC")/1e3;
%latCMITbound = ncread(froot+"/cases/PTDC_002/output/275001/MITgcm/output.nc","YC")/1e3;

g(3)=plot([-1.6e3  min(lonCMITbound(:)) min(lonCMITbound(:))],...
    [min(latCMITbound(:)) min(latCMITbound(:)) -3.4e2],...
    "color",[0.31 0.28 0.7],"linestyle",":","linewidth",4);%[246,190,0]/255
g(2)=plot([min(lonCMITbound(:)) min(lonCMITbound(:)) max(lonCMITbound(:)) max(lonCMITbound(:)) -1.6e3],...
    [-3.4e2 max(latCMITbound(:)) max(latCMITbound(:)) min(latCMITbound(:)) min(latCMITbound(:))],...
    "color",[0.31 0.28 0.7],"linestyle","-","linewidth",4);%[246,190,0]/255

%[0.4 0.4 1]
lonCMIT1=lonCMIT; latCMIT1=latCMIT;
lonCMIT1(mask==2)=NaN;
latCMIT1(mask==2)=NaN;

PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"k","linewidth",2);

BasinShape = shaperead(froot_data+"/Antarctica_DrainageBasins_IMBIE/Basins_Antarctica_v02.shp");

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

PlotLengthScale(xmin/1e3+900,ymin/1e3+230,[0 100 200 300],10);

axes;
CM2= CM;
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

%xlim([-2.5e6 -0.5e6]/1e3); ylim([-1.5e6 1.7e6]/1e3);
xlim([xmin xmax]/1e3); ylim([ymin ymax]/1e3);

box off; grid off;

% a

% % Get the color data of the object that correponds to the colorbar
% cdata = cb.Face.Texture.CData;
% % Change the 4th channel (alpha channel)
% cdata(end,:) = uint8(alphaVal*cdata(end,:));
% % Ensure that the display respects the alpha channel
% cb.Face.Texture.ColorType = "truecoloralpha";
% % Update the color data with the new transparency information
% cb.Face.Texture.CData = cdata;
% drawnow
% % Make sure that the renderer doesn"t revert your changes
% cb.Face.ColorBinding = "discrete";

% subplot("position",[0.54 0.1 0.43 0.85]); hold on;
% 
% %lima(-1400e3,-350e3,1000,"xy");
% 
% plot(lonCMIT,latCMIT,"color",[0.4 0.4 1]);
% plot(lonCMIT",latCMIT","color",[0.4 0.4 1]);
% 
% g(1) = PlotMuaMesh(CtrlVar,MUA);
% 
% plot(lonCMIT1,latCMIT1,"color",[0.8 0.6 0.8]);
% plot(lonCMIT1",latCMIT1","color",[0.8 0.6 0.8]);
% 
% g(2) = plot([min(lonCMIT_orig(:)) min(lonCMIT_orig(:)) max(lonCMIT_orig(:)) max(lonCMIT_orig(:)) min(lonCMIT_orig(:))],...
%     [min(latCMIT_orig(:)) max(latCMIT_orig(:)) max(latCMIT_orig(:)) min(latCMIT_orig(:)) min(latCMIT_orig(:))],"color",[0.4 0.4 1],"linewidth",1.5);
% 
% PlotGroundingLines(CtrlVar,MUA,GF,[],[],[],"k");
% 
% for ii=1:length(BasinShape)
%      X = double(BasinShape(ii).X)/1e3;
%      Y = double(BasinShape(ii).Y)/1e3;
%     plot(X,Y,"color",[0.4 0.4 0.4]); 
% end
% 
% %% ice fronts
% load("HeatVolumeTransport_IceFront_below400m_PTDC_001.mat");
% basins = {"PIG","TW","CR","DT"};
% color = [[51 34 136];[17 119 51];[221 204 119];[136 34 85]]/255;
% for gg=1:numel(basins)
%     basin = basins{gg};
%     switch basin
%         case "PIG"
%             xpoly = [-1699e3 -1699e3 -1530e3 -1530e3];
%             ypoly = [-380e3 -220e3 -220e3 -380e3];
%         case "TW"
%             xpoly = [-1620e3 -1620e3 -1500e3 -1500e3];
%             ypoly = [-520e3 -380e3 -380e3 -520e3];
%         case "CR"
%             xpoly = [-1580 -1485 -1450 -1450 -1610]*1e3;
%             ypoly = [-579 -657 -657 -520 -520]*1e3; 
%         case "DT"
%             xpoly = [-1610 -1485 -1450 -1450 -1625]*1e3;
%             ypoly = [-630 -630 -630 -700 -700]*1e3; 
%     end
%     addpath([froot,"/Matlab/Matlab_Functions/poly_stuff"]);
%     I = find(inpoly([section(1).xmid(:) section(1).ymid(:)],[xpoly(:) ypoly(:)]));
%     plot(section(1).xmid(I)/1e3,section(1).ymid(I)/1e3,"-","color",color(gg,:),"linewidth",4);
% end
% 
% title("");
% axis tight;
% 
% %xlim([-1.72e6 -1.4e6]/1e3); ylim([-0.72e6 -0.2e6]/1e3);
% xlim([min(MUA.coordinates(:,1)) max(MUA.coordinates(:,1))]); 
% ylim([min(MUA.coordinates(:,2)) max(MUA.coordinates(:,2))]); 
% grid on; box on;

legend(g(:),{"$\mbox{{\'U}a ice boundary}$","$\mbox{MITgcm ocean boundary}$","$\mbox{MITgcm open boundary}$"},"location","northeast","interpreter","latex","fontsize",14);

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/Figure1";
print(H,fname,"-dpng","-r400");

end