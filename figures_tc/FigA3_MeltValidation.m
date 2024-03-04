function FigA3_MeltValidation

section = ["PIG", "TW", "CR", "DT"];
sectiontitle = ["Pine Island","Thwaites","Crosson","Dotson"];
runID = ["PTDC_000","PTDC_001"];
doplots = 1;
%timestep = [0.25 0.25];
start = 1;

froot_UaMitgcm = getenv("froot_uamitgcm");

load(froot_UaMitgcm+"/Ua_InputData/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat","FB")
load(froot_UaMitgcm+"/UaMITgcm_source/example/PTDC_666/ua_custom/BoundaryCoordinates.mat");
addpath(getenv("froot_tools"));

%%%%%%%%%%%%%%%%%%%%%
%% READ MELT RATES %%
%%%%%%%%%%%%%%%%%%%%%
for ii=1:numel(runID)

    frootm = froot_UaMitgcm+"/cases/"+runID(ii)+"/output/";
    subd=dir(frootm);
    isub = [subd(:).isdir]; %# returns logical vector
    nameFolds = {subd(isub).name}';
    nameFolds(ismember(nameFolds,[".",".."])) = [];

    i1=1; i2=1; i3=1; i4=1; starttime=0; nn=0; data(ii).melt=[];
    starttime = 0;
    %starttime = 12*100;
    %nnmax = min(length(nameFolds),23*12/(12*timestep(ii)));
    
    while starttime<=datenum("20150101","yyyymmdd") && nn<numel(nameFolds)
    %while starttime<=datenum("21150101","yyyymmdd") && nn<numel(nameFolds)

        nn=nn+1;

        MITpath = frootm+nameFolds{nn}+"/MITgcm";
        MITfile = MITpath+"/output.nc";
        nstr = strlength(MITpath);
        startdate = extractBetween(MITpath,nstr-12,nstr-7);
        starttime = datenum(startdate+"01","yyyymmdd");
        T = double(ncread(MITfile,"time"));
        timestep = numel(T)/12;
        % read time epoch
        attvalue=ncreadatt(MITfile,"time","units");
        if strfind(attvalue,"seconds")
            epoch = erase(attvalue,"seconds since ");
            epochnum = datenum(epoch);
            MITTime = epochnum + T/(24*60*60);
        elseif strfind(attvalue,"days")
            epoch = erase(attvalue,"days since ");
            epochnum = datenum(epoch);
            MITTime = epochnum + T;
        else
            error("I do not recognise the time format in output.nc");
        end
        
        if starttime>datenum("31121996","ddmmyyyy")
        %if starttime>datenum("31122096","ddmmyyyy") 

            for tt = 1:numel(T)
                
                   [LON,LAT,Melt] = PlotMeltRates(MITpath,tt,[1 0 0 0 0 0 0 0]);
                   dx = LON(2,1) - LON(1,1); dy = LAT(1,2) - LAT(1,1);
                   [nx,ny] = size(Melt);
                   
                   data(ii).melt(i1,:) = Melt(:)*365.25*24*60*60/1e3; %kg/s/m2 to m/yr
                   
                   Melt = Melt*365.25*24*60*60*dx*dy/1e12; %kg/s/m2 to Gt/yr
                   Imelt = find(Melt~=0);
    
    
                   for ss=1:numel(section)
                        switch section{ss}
                            case "PIG"
                                xmin = -1699e3; xmax = -1530e3; ymin = -380e3; ymax = -220e3;
                                IPIG = find(LON(Imelt)>xmin & LON(Imelt)<xmax & LAT(Imelt)>ymin & LAT(Imelt)<ymax);
                                data(ii).section(ss).melt_int(i1) = sum(Melt(Imelt(IPIG)),"all","omitnan");
                                data(ii).section(ss).melt_mean(i1) = mean(data(ii).melt(i1,Imelt(IPIG)),2);
                                data(ii).section(ss).time(i1) = MITTime(tt);
                                i1 = i1+1;
                            case "TW"
                                xmin = -1620e3; xmax = -1500e3; ymin = -520e3; ymax = -380e3;
                                ITW = find(LON(Imelt)>xmin & LON(Imelt)<xmax & LAT(Imelt)>ymin & LAT(Imelt)<ymax);
                                data(ii).section(ss).melt_int(i2) = sum(Melt(Imelt(ITW)),"all","omitnan");
                                data(ii).section(ss).melt_mean(i2) = mean(data(ii).melt(i2,Imelt(ITW)),2);
                                data(ii).section(ss).time(i2) = MITTime(tt);
                                i2 = i2+1;
                            case "CR"
                                xpoly = [-1610 -1485 -1450 -1450 -1610]*1e3;
                                ypoly = [-580 -657 -657 -520 -520]*1e3; 
                                ICR = find(inpoly([LON(Imelt(:)),LAT(Imelt(:))],[xpoly(:) ypoly(:)]));
                                data(ii).section(ss).melt_int(i3) = sum(Melt(Imelt(ICR)),"all","omitnan");
                                data(ii).section(ss).melt_mean(i3) = mean(data(ii).melt(i3,Imelt(ICR)),2);
                                data(ii).section(ss).time(i3) = MITTime(tt);
                                i3 = i3+1;
                            case "DT"
                                xpoly = [-1610 -1485 -1450 -1450 -1625]*1e3;
                                ypoly = [-580 -657 -657 -700 -700]*1e3; 
                                IDT = find(inpoly([LON(Imelt(:)),LAT(Imelt(:))],[xpoly(:) ypoly(:)]));
                                data(ii).section(ss).melt_int(i4) = sum(Melt(Imelt(IDT)),"all","omitnan");
                                data(ii).section(ss).melt_mean(i4) = mean(data(ii).melt(i4,Imelt(IDT)),2);
                                data(ii).section(ss).time(i4) = MITTime(tt);
                                i4 = i4+1;
                        end
                   end
            end
            
            fprintf("%s: Done %i \n",runID{ii},nn);
        end
        
    end

    for nn=1:size(data(ii).melt,2)

        data(ii).meanmelt(nn) = mean(data(ii).melt(:,nn),"all","omitnan");

    end

    data(ii).meanmelt = reshape(data(ii).meanmelt,nx,ny);

end

%%%%%%%%%%%%%%%%%%
%% DATA KAITLIN %%
%%%%%%%%%%%%%%%%%%
% PIG:
% Dutrieux 2014 + Heywood 2016 estimates of PIG melting (Gt/y)
% Density of freshwater (kg/m^3)
rho_fw = 1e3;
% Density of ice (kg/m^3)
rho_ice = 917;
dutrieux_sw = [51.3, 79.7, 75.2, 37.3]*rho_ice/rho_fw;
dutrieux_mw = [49.1, 79.4, 69.2, 34.7]*rho_ice/rho_fw;
dutrieux_melt = 0.5*(dutrieux_sw+dutrieux_mw);
dutrieux_err_fac = 0.1;
heywood_melt = [40]*rho_ice/rho_fw;
heywood_err_fac = 0.4;
pig_melt_years = [1994, 2009, 2010, 2012, 2014];
pig_melt_years = datenum(num2str(pig_melt_years'),'yyyy');
pig_melt = [dutrieux_melt, heywood_melt];
pig_err = [dutrieux_melt.*dutrieux_err_fac, heywood_melt.*heywood_err_fac];

% Dotson:
% Jenkins 2018 estimates of Dotson melting (Gt/y)
dotson_melt_years = [2000, 2006, 2007, 2009, 2011, 2012, 2014, 2016];
dotson_melt_years = datenum(num2str(dotson_melt_years'),'yyyy');
dotson_melt = [25.0, 55.7, 44.4, 91.6, 53.8, 20.3, 20.9, 19.5];
dotson_err = [9.1, 15.3, 35.2, 31.6, 12.3, 6.1, 8.0, 15.4];

% PAS simulations
pig_massloss = ncread('pig_dotson_massloss.nc','pig_massloss');
dotson_massloss =  ncread('pig_dotson_massloss.nc','dotson_massloss');
time_massloss = ncread('pig_dotson_massloss.nc','time');
attvalue=ncreadatt('pig_dotson_massloss.nc',"time","units");
epoch = erase(attvalue,"seconds since ");
epochnum = datenum(epoch);
time_massloss = double(epochnum+time_massloss/(24*60*60));

%%%%%%%%%%%%%%
%% PLOTTING %%
%%%%%%%%%%%%%%
%CM = othercolor("Set14",numel(runID));

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
CM = flipdim(CM,1);

H=fig("units","inches","width",80*12/72.27,"height",45*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(2,2,"TileSpacing","tight"); n=1;
for i = [1 3]
    ax_fig(n) = nexttile(tlo_fig,i); hold on;
    n=n+1;
end
ax_fig(n) = nexttile(tlo_fig,2,[2,1]); hold on;

ss=1;
for n=[1 4]

    if contains(section{n},"PIG")
        h(3) = plot(ax_fig(ss),time_massloss,pig_massloss,'--','color',[0 0.45 0.74],'linewidth',1);       
    elseif contains(section{n},"DT")
        plot(ax_fig(ss),time_massloss,dotson_massloss,'--','color',[0 0.45 0.74],'linewidth',1); 
    end

    h(2) = plot(ax_fig(ss),data(1).section(n).time,-data(1).section(n).melt_int,"--","color",[0.3 0.3 0.3],"linewidth",1);
    h(1) = plot(ax_fig(ss),data(2).section(n).time,-data(2).section(n).melt_int,"-","color","k","linewidth",1.5);

    if contains(section{n},"PIG")
        h(4) = errorbar(ax_fig(ss),pig_melt_years,pig_melt,pig_err,'.r','LineWidth',2,'markersize',0.1);% yobs,fluxobs,"dk","markersize",10,"markerfacecolor","m"
    elseif contains(section{n},"DT")
        errorbar(ax_fig(ss),dotson_melt_years,dotson_melt,dotson_err,'.r','linewidth',2,'markersize',0.1); 
    end

    ylim(ax_fig(ss),[0 140]);
    xlim(ax_fig(ss),[datenum("01011990","ddmmyyyy") datenum("01012015","ddmmyyyy")]);
    xticks(ax_fig(ss),[datenum("01011990","ddmmyyyy") datenum("01011995","ddmmyyyy") ...
        datenum("01012000","ddmmyyyy") datenum("01012005","ddmmyyyy")...
         datenum("01012010","ddmmyyyy") datenum("01012015","ddmmyyyy")]);
        %datenum("01012009","ddmmyyyy") datenum("01012012","ddmmyyyy")...

    if ss==2       
        xlabel(ax_fig(ss),"Time");
        xticklabels(ax_fig(ss),{'1990','1995','2000','2005','2010','2015'}); 
        legend(ax_fig(ss),h(:),{"$\mbox{{\'U}aMITgcm (\textit{var\_melt})}$",...
            "$\mbox{MITgcm (\textit{ref\_melt})}$","$\mbox{MITgcm (Naughten et al., 2022)}$","$\mbox{Observations}$"},'Interpreter','latex',...
            'fontsize',12,'location','northwest');
            
    else
        xticklabels(ax_fig(ss),"");     
        
    end

    grid(ax_fig(ss),"on"); axis(ax_fig(ss),"on"); box(ax_fig(ss),"on");
    title(ax_fig(ss),sectiontitle{n});

    ss=ss+1;

end    

data(2).meanmelt(data(2).meanmelt==0)=NaN;

data(2).meanmelt(data(2).meanmelt<-100)=-100;
contourf(ax_fig(3),LON/1e3,LAT/1e3,-data(2).meanmelt,[0:5:100],"LineStyle","none"); 
%contourf(ax_fig(3),LON/1e3,LAT/1e3,data(1).meanmelt-data(2).meanmelt,[-100:1:100],"LineStyle","none"); 
shading(ax_fig(3),"flat");

frootm = froot_UaMitgcm+"/cases/"+runID(2);

restartfile = dir(frootm+"/output/199701/ua/Ua*.mat");
load(restartfile(1).folder+"/"+restartfile(1).name,"MUA","GF","CtrlVarInRestartFile");
plot(ax_fig(3),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k","linewidth",1);

uafile = dir(frootm+"/output/201501/ua/Ua*.mat");
if ~isempty(uafile)
    load(uafile(1).folder+"/"+uafile(1).name,"MUA","GF","CtrlVarInRestartFile");
end
CtrlVarInRestartFile.PlotGLs=0;
[xGL,yGL,~]=PlotGroundingLines(CtrlVarInRestartFile,MUA,GF);
plot(ax_fig(3),xGL/1e3,yGL/1e3,"-k","linewidth",0.5);


caxis(ax_fig(3),[0 100]);
%caxis(ax_fig(3),[-100 100]);
%CM = parula; CM = flipdim(CM,1);
%CM = [[linspace(0.75,0,64)';linspace(0,0.75,64)'] linspace(0.75,0,128)' linspace(0.75,1,128)'];
colormap(ax_fig(3),CM);

cb=colorbar(ax_fig(3),"Location","eastoutside","Ticks",[0:20:100]);
cb.XColor="k";
cb.YColor="k";
cb.TickLength=0.025;
cb.FontSize=16;
cb.Label.String = "$\mbox{Average basal melt 1997-2014 [m yr}^{-1}\mbox{]}$";
cb.Label.Interpreter = "latex";

% basins = DefineBasins;
% toplot(5).I = [2:44];
% toplot(4).I = [2:91];
% toplot(3).I = [48:168 1:47];
% toplot(2).I = [189:222 1:188];
% 
% for n=2:5
%     x=basins(n).X; y=basins(n).Y;
%     %I = find(inpoly([x(:) y(:)],MUA.coordinates));
%     %plot(ax_fig(5),x(I)/1e3,y(I)/1e3,"-k");
%     plot(ax_fig(5),x(toplot(n).I)/1e3,y(toplot(n).I)/1e3,"-k");
% end

axis(ax_fig(3),"equal");
xlim(ax_fig(3),[-1700 -1470]); ylim(ax_fig(3),[-700 -250]);
xlabel(ax_fig(3),"psx [km]"); ylabel(ax_fig(3),"psy [km]");

grid(ax_fig(3),"on"); axis(ax_fig(3),"on"); box(ax_fig(3),"on");
title(ax_fig(3),"$\mbox{{\'U}aMITgcm (\textit{var\_melt})}$","interpreter","latex");

ylabel(tlo_fig,"$\mbox{Basal mass loss [Gt yr}^{-1}\mbox{]}$","interpreter","latex","fontsize",16);

% [X,Y,w_b_Gourmelen,w_b_Adus] = Plot_observations_basalmelt(0);
% w_b_Gourmelen(isnan(w_b_Gourmelen))=0;
% w_b_Adus(isnan(w_b_Adus))=0;
% w_b = w_b_Gourmelen + (w_b_Gourmelen==0).*w_b_Adus;
% w_b(w_b==0)=NaN;
% 
% N=1;
% pcolor(ax_fig(6),X(1:N:end,1:N:end)/1e3,Y(1:N:end,1:N:end)/1e3,w_b(1:N:end,1:N:end)); shading(ax_fig(6),"flat");
% colormap(ax_fig(6),CM); caxis(ax_fig(6),[0 100]); 
% 
% plot(ax_fig(6),MUA.Boundary.x/1e3,MUA.Boundary.y/1e3,"-k","linewidth",1);
% plot(ax_fig(6),xGL/1e3,yGL/1e3,"-k","linewidth",0.5);
% 
% axis(ax_fig(6),"equal");
% xlim(ax_fig(6),[-1700 -1470]); ylim(ax_fig(6),[-700 -250]);
% xlabel(ax_fig(6),"psx [km]"); yticklabels(ax_fig(6),"");
% 
% grid(ax_fig(6),"on"); axis(ax_fig(6),"on"); box(ax_fig(6),"on");
% title(ax_fig(6),"$\mbox{N. Gourmelen, 2010-2016}$","interpreter","latex");
% 
% 
% legend(ax_fig(1),[h1;h2],{"$\mbox{{\'U}a-MITgcm}$";fluxcitation1},"Location","northwest","interpreter","latex");
% legend(ax_fig(4),[h3],{fluxcitation2},"Location","northwest","interpreter","latex");
% 
% datetick(ax_fig(4),"x","yyyy","keeplimits","keepticks");
% ylabel(tlo_fig,"Freshwater flux [Gt yr^{-1}]","fontsize",13);
% 
% Pos = cb.Position; cb.Position=[Pos(1)+0.015 Pos(2) Pos(3) Pos(4)];

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/FigA2_Meltrates_UaMITgcm_vs_Obs";
print(H,fname,"-dpng","-r400");

return
H=fig("units","inches","width",50*12/72.27,"height",50*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(2,2,"TileSpacing","tight");
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

MeanMelt = 0;

for ss=1:numel(section)
        
        I = find(data(1).section(ss).time(:)==datenum("01112002","ddmmyyyy"));
        Mw = -data(1).section(ss).melt_mean(I);
        
        histogram(ax_fig(ss),-data(1).section(ss).melt_mean(:),[0:1:20],"normalization","probability");
        xlim(ax_fig(ss),[0 21]); ylim(ax_fig(ss),[0 0.25]);
        grid(ax_fig(ss),"on"); box(ax_fig(ss),"on");
        title(ax_fig(ss),section{ss});

%         nmonths =[1:12*(2015-1997)];
%         CM = parula(numel(nmonths));
%         
%         for nn=nmonths
%             mm = mod(nn,12); if mm==0; mm=12; end
%             yyyy = 1997+floor(nn/13);
%             [~,I] = min(abs(data(ii).section(ss).time(:)-datenum(["01",sprintf("%02i",mm),num2str(yyyy)],"ddmmyyyy")+30));
%             M = -data(1).section(ss).melt_mean(I);
%             plot(ax_fig(ss),[M M],[0 0.25],"-","color",CM(nn,:),"linewidth",0.5);
%         end


        plot(ax_fig(ss),[Mw Mw],[0 0.25],"-r","linewidth",1.5);
        plot(ax_fig(ss),0.5*[Mw Mw],[0 0.25],"--r","linewidth",0.5);
        plot(ax_fig(ss),0.25*[Mw Mw],[0 0.25],"--r","linewidth",0.5);

        yearstoplot = [730638]; %01062000 % 735720 735781];
        CM = parula(numel(yearstoplot));
        for nn=1:numel(yearstoplot)
            I = find(data(ii).section(ss).time(:)==yearstoplot(nn));
            M = -data(1).section(ss).melt_mean(I);
            plot(ax_fig(ss),[M M],[0 0.25],"-","color",CM(nn,:),"linewidth",0.5);
        end

        %MeanMelt = MeanMelt - data(1).section(ss).melt_mean(:);
        %Time = data(1).section(ss).time(:);
        
        if ss<3
            xticklabels(ax_fig(ss),{""});
        end
        if mod(ss,2)==0
            yticklabels(ax_fig(ss),{""});
        else
            yticks(ax_fig(ss),[0:0.05:0.25]);
            yticklabels(ax_fig(ss),{"0","5","10","15","20","25"});
        end

end

xlabel(tlo_fig,"Mean melt rate [m/yr]");
ylabel(tlo_fig,"Frequency [%]");


pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
folder = ["../Figures/",runID{1}];
fname = [folder,"/MeltFrequency_",runID{1}];
print(H,fname,"-dpng","-r400");

