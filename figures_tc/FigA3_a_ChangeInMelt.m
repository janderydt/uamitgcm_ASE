function FigA3_a_ChangeInMelt

plotoption ="melt";%icethickness_speed, melt_thickness

runID="ASE_himelt";
basins = ["PIG","TW","CR","DT"];
glaciername = ["Pine Island","Thwaites","Crosson","Dotson"];


timestep1 = 2*12; % after 5-year relaxation
timestep2 = 200*12;


frootm = getenv("froot_uamitgcm")+"/cases/"+runID;
%frootm = ["/media/janryd69/mainJDeRydt/UaMITgcm_v2/cases/",runID];

if timestep1>=timestep2
   error("check timesteps"); 
end

subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = string({subd(isub).name});
nameFolds(ismember(nameFolds,[".",".."])) = [];
nfolders = numel(nameFolds);
	
output = []; kk=1; jj=1;
while kk<=timestep2 && jj<nfolders
    disp(frootm+"/output/"+nameFolds{jj});
    output(jj).filelist=dir(frootm+"/output/"+nameFolds(jj)+"/Ua/UaDefault*.mat");
    % only use 1st file
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

nstr = strlength(files{timestep1+1}.name);
year0 = extractBetween(files{timestep1+1}.name,nstr-14,nstr-11);   
month0 = extractBetween(files{timestep1+1}.name,nstr-10,nstr-9);   
if timestep1==1
    CtrlVar = CtrlVarInRestartFile;
    UserVar = UserVarInRestartFile;
    rho = F.rho; s = F.s; h = F.h; ub = F.ub; vb = F.vb; ab=F.ab;
    b = F.b; B=F.B;
    load(files{timestep1+1}.folder+"/"+files{timestep1+1}.name,"time");
end
time0 = datestr(datenum("01/"+month0+"/"+year0,"dd/mm/yyyy"),"dd/mm/yyyy");

%load(getenv("froot_tools")+"CavityStreamFunctions_"+runID+"_monthly_v2.mat");
load(getenv("froot_tools")+"HeatVolumeTransport_moving400mdraft_"+runID+".mat");
[~,I] = min(abs(datenum("01/"+month0+"/"+year0,"dd/mm/yyyy")-MITTime));
bsf.contours_20170101.xmid = section(I).xmid;
bsf.contours_20170101.ymid = section(I).ymid;
    
MUA0 = MUA; GF0=GF;
CtrlVar.PlotXYscale = 1e3;
CtrlVar.PlotGLs = 0;
[xGL1,yGL1,~] = PlotGroundingLines(CtrlVar,MUA,GF);
%Fb1 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),b);
Fab1 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ab);

%if timestep2 > length(files)
%    timestep2 = length(files);
%end

nstr = strlength(files{timestep2}.name);
year1 = string(extractBetween(files{timestep2}.name,nstr-14,nstr-11));  
month1 = string(extractBetween(files{timestep2}.name,nstr-10,nstr-9));
days1 = double(string(extractBetween(files{timestep2}.name,nstr-7,nstr-4)));
datenum("01/"+month1+"/"+year1);
time1 = datestr(datenum("01/"+month1+"/"+year1,"dd/mm/yyyy")+days1,"dd/mm/yyyy");
load(files{timestep2}.folder+"/"+files{timestep2}.name,"MUA","GF","ab");

CtrlVar.PlotGLs = 0;
[xGL2,yGL2,~] = PlotGroundingLines(CtrlVar,MUA,GF);
Fab2 = scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),ab);

[~,I] = min(abs(datenum("01/"+month1+"/"+year1,"dd/mm/yyyy")-MITTime));
bsf.contours_22141201.xmid = section(I).xmid;
bsf.contours_22141201.ymid = section(I).ymid;

xmin2 = min(MUA.coordinates(:,1)); xmax2 = max(MUA.coordinates(:,1));
ymin2 = min(MUA.coordinates(:,2)); ymax2 = max(MUA.coordinates(:,2));

[Xm,Ym] = ndgrid([xmin2:500:xmax2],[ymin2:500:ymax2]);

S_insar = shaperead(getenv("froot_data")+"/GroundingLine/MEaSUREs_INSAR/InSAR_GL_Antarctica_v02.shp");
S_insar_mililloTW = shaperead(getenv("froot_data")+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_Thwaites/Thwaites_2016-2017.shp");
S_insar_mililloDC = shaperead(getenv("froot_data")+"/PIG_Thwaites_DotsonCrosson/Milillo_GL_CrossonDotson/Grounding_Lines_CSK_Pope_Smith_Kohler.shp");

%% PLOTTING

H=fig('units','inches','height',70*12/72.27,'width',60*12/72.27,'fontsize',16,'font','Helvetica','visible','off');

subplot("position",[0 0 1 1]); hold on;

melt1 = Fab1(Xm,Ym);
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y]));
melt1(I) = NaN; melt1(melt1==0)=NaN; melt1(melt1<-200)=-200;
contourf(Xm,Ym,-melt1,[0:10:200],"LineStyle","none"); 
%CM1=othercolor('RdYlBu_11b',4*32+30); CM1 = flipdim(CM1,1);
%CM2=othercolor('RdYlBu_11b',16*32); CM2 = flipdim(CM2,1);
%CM=[CM1(1:65,:); CM2(7*32:end,:)]; CM= flipdim(CM,1);
%CM = othercolor('RdYlBu_11b',29);
CM = [0.4941    0.1843    0.5569;...
    0.5496    0.1936    0.5173;...
    0.6051    0.2029    0.4777;...
    0.6606    0.2123    0.4381;...
    0.7161    0.2216    0.3985;...
    0.7618    0.2315    0.3680;...
    0.7919    0.2424    0.3522;...
    0.8201    0.2534    0.3381;...
    0.8120    0.2503    0.3421;...
    0.8344    0.2684    0.3393;...
    0.9414    0.3656    0.3319;...
    0.9710    0.4627    0.3735;...
    0.9867    0.5547    0.4087;...
    0.9946    0.6400    0.4428;...
    0.9983    0.7197    0.4835;...
    0.9996    0.7900    0.5290;...
    0.9999    0.8634    0.5882;...
    0.9989    0.9432    0.6700;...
    0.9925    0.9683    0.7240;...
    0.9797    0.9854    0.7838;...
    0.9590    0.9954    0.8464;...
    0.9267    0.9982    0.9073;...
    0.8848    0.9970    0.9506;...
    0.8119    0.9864    0.9796;...
    0.7008    0.9600    0.9952;...
    0.5297    0.8851    0.9990;...
    0.3731    0.7412    1.0000;...
    0.2433    0.5450    1.0000];
CM = [interp1(1:size(CM,1),CM(:,1),linspace(1,size(CM,1),21))',...
    interp1(1:size(CM,1),CM(:,2),linspace(1,size(CM,1),21))',...
    interp1(1:size(CM,1),CM(:,3),linspace(1,size(CM,1),21))'];
%...
%    0.1500    0.3000    1.0000];
colormap(flipdim(CM,1));
cb=colorbar("Position",[0.12 0.3 0.017 0.35],"Location","east","AxisLocation","out","Ticks",[0:25:200],...
    "TickLabels",["0" "" "50" "" "100" "" "150" "" ">200"]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.04;
cb.FontSize=14;
cb.Label.String = "Basal melt rate [m/a]";
caxis([0 200]);
hold on;
%plot(xB.bsf_20170101,yB.bsf_20170101,'-k');
plot(xGL1,yGL1,'-','color',[0.5 0.5 0.5],'linewidth',1);
plot(MUA.Boundary.x,MUA.Boundary.y,'-','color',[0.5 0.5 0.5],'linewidth',1)
plot(bsf.contours_20170101.xmid,bsf.contours_20170101.ymid,'ok','markersize',2,'MarkerFaceColor',[0.2 0.2 0.2]);
axis equal; 
xticklabels(gca,{}); yticklabels(gca,{});
set(gca,"Visible","off");

%% depth average flow
% quiver(X,Y,u_bt.*xymask,v_bt.*xymask,"color","k","AutoScaleFactor",10);

%subplot("position",[0.251 0 0.4 1]); hold on;
melt2 = Fab2(Xm,Ym); melt2(melt2==0)=NaN;
I = find(~inpoly([Xm(:) Ym(:)],[MUA.Boundary.x MUA.Boundary.y]));
melt2(I) = NaN; melt2(melt2<-200)=-200;
contourf(Xm+0.16e6,Ym,-melt2,[0:10:200],"LineStyle","none"); 
caxis([0 200]);
hold on
%plot(xB.bsf_20170101+0.18e6,yB.bsf_20170101,'-k');
%plot(xGL.bsf_20170101,yGL.bsf_20170101,'-k');
plot(xGL2+0.16e6,yGL2,'-','color',[0.5 0.5 0.5],'linewidth',1);
plot(MUA.Boundary.x+0.16e6,MUA.Boundary.y,'-','color',[0.5 0.5 0.5],'linewidth',1)
plot(bsf.contours_22141201.xmid+0.16e6,bsf.contours_22141201.ymid,'ok','markersize',2,'MarkerFaceColor',[0.2 0.2 0.2]);
ylim([-7e5 -2.25e5]); xlim([-1.69e6 -1.45e6+0.2e6]);

pos = get(H,'Position');
set(H,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
fname = ['../Figures/FigS3_1_Delta_melt_2D'];
print(H,fname,'-dpng','-r400');

return

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

    switch plotoption
        
    
        case "melt"
            
            PlotNodalBasedQuantities_JDR(ax_fig(bb),MUA.connectivity,MUA.coordinates/1e3,deltawct); hold on;
    
            %[M,c]=contour(ax_fig(bb),Xm/1e3,Ym/1e3,wctm,[0:200:1000]); c.Color=[0.7 0.7 0.7];
            %clabel(M,c,[0:200:1000],'color',[0.5 0.5 0.5],'labelspacing',500);
    
            Bcont(Bcont<-1600)=-1600; Bcont(Bcont>-200)=-200;
            scale = (-Bcont+max(Bcont,[],'omitnan'))/(max(Bcont,[],'omitnan')-min(Bcont,[],'omitnan'));
            Mask = GF.node;
            alphaMask = GF.node; I = find(alphaMask<1);
            alphaMask(I)=0;
    
            patch(ax_fig(bb),'faces',MUA.connectivity,'vertices',MUA.coordinates/1e3,...
            'FaceVertexCData',Mask,'CDataMapping','scaled','EdgeColor','none','FaceColor',[0.5 0.5 0.8],...
            'FaceAlpha','flat','FaceVertexAlphaData',alphaMask.*scale);
    
            [M2,c2]=contour(ax_fig(bb),Xm/1e3,Ym/1e3,Bm,[-1600:100:-600]); c2.Color=[0.5 0.5 1];
            clabel(M2,c2,[-1600:200:-200],'color',[0.3 0.3 1],'labelspacing',500); 
     
            
            for ii=1:length(S_insar)
                date=S_insar(ii).DATE1;
                year_insar = date(1:4);
                [x_insar,y_insar] = ll2psxy(S_insar(ii).Y,S_insar(ii).X,-71,0);
                if str2double(year_insar)>=2011
                    ginsar=plot(ax_fig(bb),x_insar/1e3,y_insar/1e3,'color',[255, 165, 0]/255,'linewidth',2.5); 
                end
            end

            for ii=1:length(S_insar_mililloTW)
                plot(ax_fig(bb),S_insar_mililloTW(ii).X/1e3,S_insar_mililloTW(ii).Y/1e3,'color',[255, 165, 0]/255,'linewidth',2.5);  
            end

            for ii=1:length(S_insar_mililloDC)
                plot(ax_fig(bb),S_insar_mililloDC(ii).X/1e3,S_insar_mililloDC(ii).Y/1e3,'color',[255, 165, 0]/255,'linewidth',2.5);                
            end
            
            plot(ax_fig(bb),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,'-k');
            
            CtrlVar.PlotGLs=0; CtrlVar.PlotXYscale = 1e3;
            [xGL,yGL,~]=PlotGroundingLines(CtrlVar,MUA,GF);
            gl_200 = plot(ax_fig(bb),xGL/CtrlVar.PlotXYscale,yGL/CtrlVar.PlotXYscale,'m','linewidth',1.5) ;
            [xGL0,yGL0,~]=PlotGroundingLines(CtrlVar,MUA0,GF0);
            gl_0 = plot(ax_fig(bb),xGL0/CtrlVar.PlotXYscale,yGL0/CtrlVar.PlotXYscale,'k','linewidth',1.5) ;

            if bb==1
                legend(ax_fig(bb),[ginsar,gl_0,gl_200],{'$\mbox{InSAR 2011-2017}$','$\mbox{{\''U}a-MITgcm year 0}$','$\mbox{{\''U}a-MITgcm year 200}$'},'location','northeast','fontsize',14,'interpreter','latex');
            end

            axis(ax_fig(bb),'equal'); box(ax_fig(bb),'on');
            xlim(ax_fig(bb),[xmin xmax]); ylim(ax_fig(bb),[ymin ymax]);
    
            %set([ax1],'Position',[0.1 0.1 0.85 0.85]);
    
            %% Give each one its own colormap
            CM = othercolor('RdYlBu11',23); CM = CM(12:23,:); CM(2,:)=[1 1 1];
            colormap(ax_fig(bb),CM); caxis(ax_fig(bb),[-100 750]);
    
            if contains(basins{bb},'PIG')
                cb1=colorbar(ax_fig(bb),'Position',[0.08 0.8 0.10 0.015],'Location','northoutside','AxisLocation','in','Ticks',[-100 0:250:750]);
                cb1.XColor='k';
                cb1.YColor='k';
                cb1.TickLength=0.04;
                cb1.FontSize=14;
                cb1.Label.String = '\DeltaH [m]';
            end
    
            %CM2 = othercolor('Bu_10',6); 
            %colormap(ax2,CM2); caxis(ax2,[-1600 -600]);
            
    %            title('UaMITgcm');
            if bb==1
                    ylabel(ax_fig(bb),'psy [km]','Fontsize',14);
            elseif bb==2
                    xlabel(ax_fig(bb),'psx [km]','Fontsize',14);
            end

            title(ax_fig(bb),basinname);
            
        otherwise 
            error("case unknown")
    
    end
            
end
    
pos = get(Hfig,"Position");
set(Hfig,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);

folder = ["../Figures/",runID];
if ~exist(folder)
    mkdir(folder);
end
fname = [folder,"/Fig4_ChangeInWCT"];
print(Hfig,fname,"-dpng","-r400");