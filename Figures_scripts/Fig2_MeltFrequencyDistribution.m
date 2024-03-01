function Fig2_MeltFrequencyDistribution

%% Plotting routine for Fig 2 in https://doi.org/10.5194/egusphere-2023-1587
%% Requires the uamitgcm toolbox (https://github.com/janderydt/uamitgcm_tools)
%% and Ua Utilities (https://github.com/GHilmarG/UaSource/tree/beta/UaUtilities)

% Initialize UaMitgcm case directory 
froot_data = getenv("froot_uamitgcm");

% Load UaMITgcm toolbox
addpath(getenv("froot_tools"));

% Basins to plot
basins = ["PIG", "TW", "CRDT"];
basintitle = ["Pine Island","Thwaites","Crosson & Dotson"];

%% Gather data
% Load geometry data
load(froot_data+"/Ua_InputData/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat","FB");
load(froot_data+"/UaMITgcm_source/example/PTDC_666/ua_custom/BoundaryCoordinates.mat");

% Read melt rates for experiments in runID
runID = ["PTDC_001","PTDC_002_v1","PTDC_003"];

for id = 1:numel(runID)
    % HeatVolumeTransport files are produced with the GenerateTransportFiles.m function in the
    % UaMITgcm toolbox
    load("heatvolumetransport_icefront_below400m_"+runID(id)+".mat");
    if contains(runID(id),["PTDC_001","PTDC_000"])
        I = find(datenum("19970101","yyyymmdd")<=MITTime & MITTime<datenum("20150101","yyyymmdd"));
    else
        I = find(datenum("20160101","yyyymmdd")<=MITTime & MITTime<datenum("20170101","yyyymmdd"));
    end
    %PIG
    data(id).basins(1).melt = integral2D.PIG.monthly.melt_integral(I)./integral2D.PIG.monthly.ISarea_integral(I)*1e9;
    %TW
    data(id).basins(2).melt = integral2D.TW.monthly.melt_integral(I)./integral2D.TW.monthly.ISarea_integral(I)*1e9;
    %DC
    data(id).basins(3).melt = (integral2D.DT.monthly.melt_integral(I)+integral2D.CR.monthly.melt_integral(I))./...
        (integral2D.DT.monthly.ISarea_integral(I)+integral2D.CR.monthly.ISarea_integral(I))*1e9;
    data(id).time = MITTime(I);
end

%% Plotting
H=fig("units","inches","width",100*12/72.27,"height",30*12/72.27,"fontsize",14,"font","Helvetica");

tlo_fig = tiledlayout(1,3,"TileSpacing","compact");
for i = 1:3
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

style = ["-","-","--"];

for id=1:numel(runID)
        
    for ss=1:numel(basins)

        if contains(runID(id),["PTDC_001","PTDC_000"])
            g(id)=histogram(ax_fig(ss),data(id).basins(ss).melt(:),[0:1:20],"normalization","probability");
            M = mean(data(id).basins(ss).melt(:));
            h=plot(ax_fig(ss),[M M],[0 0.5],"-m","linewidth",1);
        else
            M = mean(data(id).basins(ss).melt(:));
            g(id)=plot(ax_fig(ss),[M M],[0 0.5],"-k","linewidth",2.5,"linestyle",style(id));
        end
        xlim(ax_fig(ss),[0 21]); ylim(ax_fig(ss),[0 0.25]);
        grid(ax_fig(ss),"on"); box(ax_fig(ss),"on");
        title(ax_fig(ss),basintitle{ss});

        yticks(ax_fig(ss),[0:0.05:0.25]);
        if ss==1
            yticklabels(ax_fig(ss),{"0","5","10","15","20","25"});
        else
            yticklabels(ax_fig(ss),{""});
        end   

    end

end

legend(ax_fig(1),[g(1) h g(2) g(3)],["1997-2014","mean 1997-2014","\it hi\_melt \rm(11/2002)","\it av\_melt \rm(01/1998)"],"location","northwest");

xlabel(tlo_fig,"Mean melt rate [m/yr]","fontsize",16);
ylabel(tlo_fig,"Frequency [%]","fontsize",16);

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
folder = ["../Figures/"];
fname = [folder+"/Fig2_MeltFrequency"];
print(H,fname,"-dpng","-r400");

