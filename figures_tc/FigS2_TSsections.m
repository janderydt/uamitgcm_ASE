function FigS2_TSsections

variables_to_plot = "TS"; %options: TS, VEL, gradRho, gradT, gradS
sections = ["DT","PIG"];

runID = "ASE_himelt";
timestep = 1;% in number of months

rhoConst = 1024;
si2dbar = 1e-4;
gravity = 9.81;
a0 = -0.0573; a1 = 0; a2 = 0; c0 = 0.0832; b0 = -7.53e-4; % coefficients for local freezing temp of water

start = 12*5; % start time in months
dt = 12*25; % time diff between tiles [in months]
iterations = 4; % number of slices to include

%interp = dir([frootm,"/Data/GriddedInterpolants_sBh_Bedmap2.mat"]);
%load([interp.folder,"/",interp.name],"FB");
load(getenv('froot_uamitgcm')+"/Ua_InputData/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat","FB")
load(getenv('froot_uamitgcm')+"/cases/ASE_varmelt/ua_custom/BoundaryCoordinates.mat");

frootm = getenv('froot_uamitgcm')+"/cases/"+runID;
subd=dir(frootm+"/output/");
isub = [subd(:).isdir]; %# returns logical vector
nameFolds = string({subd(isub).name});
nameFolds(ismember(nameFolds,[".",".."])) = [];
	
% listing = [];
% for jj=1:length(nameFolds)
%         listing{jj}=dir([frootm,"/Ua_data/",nameFolds{jj},"/0006*.mat"]);
% end
output = []; kk=1;
for jj=1:start+dt*iterations%length(nameFolds)
    disp(frootm+"/output/"+nameFolds(jj));
    output(jj).filelist=dir(frootm+"/output/"+nameFolds(jj)+"/MITgcm/output.nc");
    for ii=1:length(output(jj).filelist)
        listing{kk}.folder = output(jj).filelist(ii).folder;
        listing{kk}.name = output(jj).filelist(ii).name;
        kk = kk+1;
    end
end

% if start+iterations>length(listing)
%     error(["Data record (",num2str(length(listing)),") too short for start at ",...
%         num2str(start)," and ",num2str(iterations)," iterations"]);
% end

%%%%%%%%%%%%%%
%% PLOTTING %%
%%%%%%%%%%%%%%
folder = "../Figures/"+runID;
if ~exist(folder)
    mkdir(folder);
end

G=fig("units","inches","width",90*12/72.27,"height",45*12/72.27,"fontsize",12,"font","Helvetica");

tlo_fig = tiledlayout(numel(sections),iterations,'TileSpacing','compact');
for i = 1:numel(sections)*iterations
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

% simulation month
Imonth = reshape([1:numel(listing)*timestep],timestep,numel(listing));
MonthsToPlot = [start:dt:iterations*dt]

for it=1:numel(MonthsToPlot)
    [Itimestamp(it),IMITfile(it)] = find(Imonth==MonthsToPlot(it));
end

for it=1:numel(MonthsToPlot)

    Tii = Itimestamp(it);
  
	MITpath = listing{IMITfile(it)}.folder;
    MITfile = MITpath+"/"+listing{IMITfile(it)}.name;

    disp("Plotting "+MITpath+" ("+num2str(it)+" out of "+string(iterations)+")");
    
    % read time epoch
    Time = double(ncread(MITfile,"time"));
    attvalue=ncreadatt(MITfile,"time","units");
    if strfind(attvalue,"seconds")
        epoch = erase(attvalue,"seconds since ");
        epochnum = datenum(epoch);
        Time = epochnum + Time/(24*60*60);
    elseif strfind(attvalue,"days")
        epoch = erase(attvalue,"days since ");
        epochnum = datenum(epoch);
        Time = epochnum + Time;
    else
        error("I do not recognise the time format in output.nc");
    end
   
    Uafiles=dir(listing{IMITfile(it)}.folder(1:end-6)+"/Ua/UaDefault*.mat");
    if isempty(Uafiles)
        Uafiles = dir(frootm+"/ua_custom/"+runID+"*.mat");
        Uafile = Uafiles.folder+"/"+Uafiles.name;
    else
        Uafile = Uafiles(Tii).folder+"/"+Uafiles(Tii).name;
    end

    Theta_full = ncread(MITfile,"THETA");
    [In,Jn,Kn,T]=size(Theta_full);
    Theta = squeeze(Theta_full(:,:,:,Tii));
    Salt_full = ncread(MITfile,"SALT");
    Salt = squeeze(Salt_full(:,:,:,Tii));
    Salt(Salt==0)=NaN;
    fid=fopen(MITpath+"/bathymetry.shice","r","ieee-be"); 
    B=fread(fid,"real*8");fclose(fid);

	lon=double(ncread(MITfile,"XC"));
	lat=double(ncread(MITfile,"YC")); 
	z=-double(ncread(MITfile,"Z"));
    depthC = cumsum(ncread(MITfile,"drC"));
    zE(1) = 0;
    for kk=2:numel(z)
        zE(kk) = zE(kk-1)+2*(z(kk-1)-zE(kk-1));
    end
    zE(end+1)=depthC(end);
    dz = zE(2:end)-zE(1:end-1);

	[m,n]=size(lon);
	if n>1
	    lon = lon(:,1);
	    lat = lat(1,:);
    end
    dlon = lon(2)-lon(); dlat = lat(2)-lat(1);
	[LAT,LON,Z]=meshgrid(lat,lon,z);
	[LON2,LAT2] = ndgrid(lon,lat);

	load(Uafile,"MUA","b","s","B","GF");
    x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

    %lonsect = [-1700e3:0.5e3:-1580e3]; latsect = (-200e3+700e3)/(-1580e3+1700e3)*(lonsect+1700e3)-700e3;

    if it==1

        for ss=1:numel(sections)

            switch sections(ss)
                case "PIG"
                    coordinates_input = 1e6*[-1.5859   -0.24;  -1.619   -0.3686];

                case "DT"
                    coordinates_input = 1e6*[-1.475   -0.662;  -1.495   -0.67; -1.5306 -0.6529];
                    
                otherwise
                    error("Unknown case");
            end
            xnodes =[]; ynodes = [];
            for sn=1:size(coordinates_input,1)-1
                x1 = coordinates_input(sn,1); y1 = coordinates_input(sn,2);
                x2 = coordinates_input(sn+1,1); y2 = coordinates_input(sn+1,2);
                L_temp = sqrt((x2-x1).^2+(y2-y1).^2); dL_temp = 0.5e3; n = round(L_temp/dL_temp);
                xnodes = [xnodes linspace(x1,x2,n+1)]; ynodes = [ynodes linspace(y1,y2,n+1)]; 
            end
            sect(ss).lonsect = [xnodes(1:end-1)' xnodes(2:end)']; sect(ss).latsect = [ynodes(1:end-1)' ynodes(2:end)']; 
            dL = sqrt((sect(ss).lonsect(:,2)-sect(ss).lonsect(:,1)).^2+(sect(ss).latsect(:,2)-sect(ss).latsect(:,1)).^2);
            sect(ss).d = cumsum(dL);
            sect(ss).epar = [(sect(ss).lonsect(:,2)-sect(ss).lonsect(:,1))./dL,(sect(ss).latsect(:,2)-sect(ss).latsect(:,1))./dL];
            sect(ss).eperp = [-sect(ss).epar(:,2),sect(ss).epar(:,1)];      
            lonmid = 0.5*(sect(ss).lonsect(:,1)+sect(ss).lonsect(:,2)); 
            latmid = 0.5*(sect(ss).latsect(:,1)+sect(ss).latsect(:,2)); 
            sect(ss).lonmid = lonmid(:);
            sect(ss).latmid = latmid(:);
    
            I = inpoly([lonmid latmid],MeshBoundaryCoordinates);
            sect(ss).J = find(I==1);
            sect(ss).Ualonsect = lonmid(sect(ss).J)'; sect(ss).Ualatsect = latmid(sect(ss).J)';
            sect(ss).LONSECT = repmat(lonmid',numel(z),1); 
            sect(ss).LATSECT = repmat(latmid',numel(z),1); 
            sect(ss).ZSECT = repmat(z(:),1,numel(lonmid));

            d = 0;
            for jj=2:length(sect(ss).lonsect)
                d(jj) = d(jj-1) + sqrt((sect(ss).lonsect(jj)-sect(ss).lonsect(jj-1)).^2+(sect(ss).latsect(jj)-sect(ss).latsect(jj-1)).^2);
            end
            sect(ss).d2 = d;
            [sect(ss).D,sect(ss).Z] = ndgrid(d/1e3,-z);

        end

            Time0 = Time;

    end

    nstr = strlength(Uafile);
    yyyymm = extractBetween(Uafile,nstr-14,nstr-9);
    if it==1
        yyyymm0 = yyyymm;
    end

    months = between(datetime(yyyymm0+"01","inputformat","yyyyMMdd"),datetime(yyyymm+"01","inputformat","yyyyMMdd"),'years');
    m = split(months,{'years'});

    titletime = string(m);

    Fb = TriScatteredInterp(x,y,b);
    Fs = TriScatteredInterp(x,y,s);

    for ss=1:numel(sections)

        Tsect = interp3(LAT,LON,Z,Theta,sect(ss).LATSECT(:),sect(ss).LONSECT(:),sect(ss).ZSECT(:)); 
        sect(ss).Tsect = reshape(Tsect,size(sect(ss).LONSECT));
        Ssect = interp3(LAT,LON,Z,Salt,sect(ss).LATSECT(:),sect(ss).LONSECT(:),sect(ss).ZSECT(:)); 
        sect(ss).Ssect = reshape(Ssect,size(sect(ss).LONSECT));
        sect(ss).bsect = Fb(sect(ss).Ualonsect,sect(ss).Ualatsect);
        sect(ss).ssect = Fs(sect(ss).Ualonsect,sect(ss).Ualatsect);

        if it==1
           sect(ss).bsect_init = sect(ss).bsect;
           sect(ss).ssect_init = sect(ss).ssect;
           sect(ss).d_init = sect(ss).d2(sect(ss).J);
        end
        
        sect(ss).Bsect = FB(sect(ss).lonmid,sect(ss).latmid);

        %% PLOTTING
        pcolor(ax_fig((ss-1)*numel(MonthsToPlot)+it),sect(ss).D,sect(ss).Z,sect(ss).Tsect');
        shading(ax_fig((ss-1)*numel(MonthsToPlot)+it),"flat");

        if it>1
            plot(ax_fig((ss-1)*numel(MonthsToPlot)+it),[sect(ss).d_init/1e3 flipdim(sect(ss).d_init/1e3,2)],...
                [sect(ss).bsect_init flipdim(sect(ss).ssect_init,2)],"--","color",[0 0 0],'linewidth',1);
        end

        
        %Tmin = 1.5; Tmax = 4;
        if contains(sections(ss),"PIG")
            [C,h]=contour(ax_fig((ss-1)*numel(MonthsToPlot)+it),sect(ss).D,sect(ss).Z,sect(ss).Ssect',[34 34.2 34.3 34.4 34.5 34.6],"color","w","linewidth",1.5);
            clabel(C,h,[34.2 34.4 34.6],"fontsize",12,"color","w","labelspacing",1e6);
            Tmin = -0.75; Tmax = 1;
            ylim(ax_fig((ss-1)*numel(MonthsToPlot)+it),[-1.3e3 400]);
        elseif contains(sections(ss),"DT")
            [C,h]=contour(ax_fig((ss-1)*numel(MonthsToPlot)+it),sect(ss).D,sect(ss).Z,sect(ss).Ssect',[34 34.2 34.3 34.4 34.5 34.6],"color","w","linewidth",1.5);
            clabel(C,h,[34 34.2 34.3 34.4 34.5 34.6],"fontsize",12,"color","w","labelspacing",0.5e5);
            Tmin = -0.75; Tmax = 0.5;
            ylim(ax_fig((ss-1)*numel(MonthsToPlot)+it),[-1.8e3 400]);
        end

%         dCM = ceil(abs(Tmin)*20);
%         CM1 = othercolor("BuOrR_14",2*dCM); CM1 = CM1(1:dCM,:);
%         dCM = ceil(3/4*abs(Tmax)*20);
%         CM2 = othercolor("BuOrR_14",2*dCM); CM2 = CM2(dCM:end,:);
%         dCM = ceil(abs(Tmax)*20)-dCM;
%         CM3 = [linspace(CM2(end,1),152/255,dCM+1)' linspace(CM2(end,2),68/255,dCM+1)' linspace(CM2(end,3),158/255,dCM+1)'];
%         CM = [CM1;CM2;CM3(2:end,:)];
        CM = flipdim(othercolor('Spectral8',32),1);
        CM = [CM;linspace(CM(end,1),95/255,5)' linspace(CM(end,2),0,5)' linspace(CM(end,3),160/255,5)'];
        colormap(ax_fig((ss-1)*numel(MonthsToPlot)+it),CM);
        clim(ax_fig((ss-1)*numel(MonthsToPlot)+it),[Tmin Tmax]);
    
        patch(ax_fig((ss-1)*numel(MonthsToPlot)+it),[sect(ss).d(sect(ss).J)/1e3; flipdim(sect(ss).d(sect(ss).J)/1e3,1)],...
            [sect(ss).bsect flipdim(sect(ss).ssect,2)],"w"); 
    
        P=patch(ax_fig((ss-1)*numel(MonthsToPlot)+it),[sect(ss).d(1)/1e3 sect(ss).d'/1e3 sect(ss).d(end)/1e3],[-1900 sect(ss).Bsect(:)' -1900],"w");
        P.FaceColor=[0.85 0.85 0.85];

        xlim(ax_fig((ss-1)*numel(MonthsToPlot)+it),[1 max(sect(ss).d)/1e3-2]);
    
        grid(ax_fig((ss-1)*numel(MonthsToPlot)+it),"on"); 
        box(ax_fig((ss-1)*numel(MonthsToPlot)+it),"on");
        set(ax_fig((ss-1)*numel(MonthsToPlot)+it), "xdir", "reverse" );
    
        %plot(ax_fig((ss-1)*numel(MonthsToPlot)+it),[0 100],[-400 -400],"-k","linewidth",1.5);
        yticks(ax_fig((ss-1)*numel(MonthsToPlot)+it),[-1800:200:200]); 
        
        if it==1
            ylabel(ax_fig((ss-1)*numel(MonthsToPlot)+it),"depth [m]");
            yticklabels(ax_fig((ss-1)*numel(MonthsToPlot)+it),["","-1600","","-1200","","-800","","-400","","0",""]);
        else
            yticklabels(ax_fig((ss-1)*numel(MonthsToPlot)+it),[]);
        end 

        if it==numel(iterations)
             
            cb=colorbar(ax_fig((ss-1)*numel(MonthsToPlot)+it));
            cb.Location="east";
            cb.Position=[0.92 mod(ss,2)*0.435+0.12 0.01 0.36];
            cb.AxisLocation="out";
            set(cb,"Ticks",[Tmin:0.25:Tmax]);
            %set(cb,"TickLabels",{"-1.7","","-1.1","","-0.5","","0.1","","0.7","","1.3"});
            cb.XColor="k";
            cb.YColor="k";
            cb.TickLength=0.03;
            cb.FontSize=12;
            cb.Label.String = '$\mbox{Temperature [}^{\circ}\mbox{C]}$';
            cb.Label.Interpreter = 'Latex';
            cb.Label.FontSize=14;

        end
        
        xlabel(tlo_fig,"Distance [km]"); 
        
        if ss==1
            title(ax_fig((ss-1)*numel(MonthsToPlot)+it),titletime+" years");
        end

    end

end

pos = get(G,"Position");
set(G,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/FigS2_sections_"+string(runID);   
print(G,fname,"-dpng","-r400");

