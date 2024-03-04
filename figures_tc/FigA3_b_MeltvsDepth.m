function FigS3_b_MeltvsDepth


runID = "ASE_himelt";

frootm = getenv("froot_tools");

basins = ["PIG","TW","CR","DT"];
glaciernames= ["Pine Island","Thwaites","Crosson","Dotson"];

years = [2017 2215];
maxdepth=-25; mindepth=1725;
%colors = [0.93 0.69 0.13;0 0 0.9];
%colors = [133 87 35;0 0 255]/255;
%colors = [0.93 0.69 0.13;0 0.45 0.74];
colors = [0.3 0.3 0.3;0.2 0.6 0.5];

H=fig("units","inches","width",60*12/72.27,"height",40*12/72.27,"fontsize",12,"font","Helvetica");

tlo_fig = tiledlayout(2,2,'TileSpacing','compact');
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

for bb = 1:numel(basins)

    for jj=1:numel(years)

        filetoread = frootm+"Melt_vs_Draft_"+runID+"_"+basins(bb)+"_"+string(years(jj))+".mat";

        if ~exist(filetoread,"file")
            
            cd(frootm);

            Melt_vs_IceGeometry_v2(basins(bb),runID,years(jj),1/12,0);

            cd("/mnt/SSD1/Documents/Projects/2022_IceOceanCoupling_PTDC/JGR/Matlab_Functions");

        end

    end

end


for bb = 1:numel(basins)

    for jj=1:numel(years)

        filetoread = frootm+"Melt_vs_Draft_"+runID+"_"+basins(bb)+"_"+string(years(jj))+".mat";

        load(filetoread);
        
        b(jj)=bar(ax_fig(bb),bincenters,-median_melt);
        set(b(jj),'FaceAlpha',0.5);
        plot(ax_fig(bb),[-400 -400],[0 200],'--k');
        grid(ax_fig(bb),"on");
        box(ax_fig(bb),"on");
        ylim(ax_fig(bb),[0 120])
        title(ax_fig(bb),glaciernames(bb));

        if bb==1
            xticklabels(ax_fig(bb),"");
            legend(b(:),["year 0","year 200"],"Location","northwest");
        elseif bb==2
            xticklabels(ax_fig(bb),"");
            yticklabels(ax_fig(bb),"");
        elseif bb==4
            yticklabels(ax_fig(bb),"");
        end

    end

end

ylabel(tlo_fig,"Average basal met rate [m/a]");
xlabel(tlo_fig,"Ice-shelf draft [m]");

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/S3_Melt_vs_draft"+runID;   
print(H,fname,"-dpng","-r400");

%% Box plot

H=fig('units','inches','width',50*12/72.27,'height',70*12/72.27,'fontsize',14,'font','Helvetica');

tlo_fig = tiledlayout(2,2,'TileSpacing','tight');
for i = 1:4
    ax_fig(i) = nexttile(tlo_fig); hold on;
end

for bb = 1:numel(basins)

    filetoread = frootm+"geometry_melt_"+basins(bb)+"_"+runID+".mat";

    if ~exist(filetoread,"file")
        
        error("file "+filetoread+" does not exist");

    else

        load(filetoread);

    end

    times = [data.Time]; 
    hn(bb)=0;

    for jj=1:numel(years)

        [~,nn] = min(abs(times-datenum(['0101',num2str(years(jj))],'ddmmyyyy')));

        data2=[]; data3=[];
        for ii=nn:nn+12
            %for jj=1:length(J)
                %K = data(ii).DepthIntervals(J(jj)).Ind;
                data3 = [data3; -data(ii).b];
                data2 = [data2; -data(ii).ab];
            %end  
        end
    
        Lbin = 50;
        Nbins= (mindepth-maxdepth)/Lbin;
        disp(['Samples: ',num2str(length(data3(:)))]);
        data3=[data3;[-mindepth:25:-maxdepth]']; data2=[data2;0*[-mindepth:25:-maxdepth]'-10];
        
        %[N,edges,bin] = histcounts(data3(:),Nbins,'BinLimits',[-mindepth-Lbin/2 -maxdepth+Lbin/2]);
        [N,edges,bin] = histcounts(data3(:),Nbins,"BinLimits",[-mindepth -maxdepth]);
        I= find(bin==0); bin(I)=[]; data2(I)=[];
%         I= find(N<3); J=[];
%         for ii=1:numel(I)
%             J = [J; find(bin==I(ii))];
%         end
%         bin(J)=[]; data2(J)=[];
        
        %PD_b =  fitdist(data3(:),'kernel','kernel','epanechnikov');
        
        %x_values = [-1600:10:0];
        %y = pdf(PD_b,x_values);
        
        %yyaxis left 
        %plot(x_values,y,'LineWidth',2);
        
        bincenters = (edges(bin+1)+edges(bin))/2; BPlabels=[]; kk=1;
        for ii=1:Nbins
            if N(ii) > 0
               if mod(ii+3,ceil(Nbins/10))==0
                   BPlabels{kk}=num2str((edges(ii+1)+edges(ii))/2,4);
               else
                   BPlabels{kk}='';
               end
               kk=kk+1;
            end  
        end

        boxplot(ax_fig(bb),data2,bincenters,'outliersize',2,'labels',BPlabels,'symbol',''); kk=0;
        
        h = findobj(ax_fig(bb),'tag','Box');
        for j=1:length(h)-hn(bb)
            h(j).Visible='off';
            hlegend(jj)=patch(ax_fig(bb),get(h(j),'XData')+[-0.15 -0.15 0.15 0.15 -0.15],get(h(j),'YData'),colors(jj,:),'EdgeColor',[0.3 0.3 0.3],'FaceAlpha',.3);
        end 
        
        legendlabel{jj}=num2str(years(jj));
        
        h = findobj(ax_fig(bb),'tag','Median');
        set(h,'linestyle','none');
        for ii=1:length(h)-hn(bb)
            BinCenter(ii) = (h(ii).XData(2)+[h(ii).XData(1)])/2;
            Median(ii) = h(ii).YData(1);
        end
        
        N = flipdim(N(:),1); I=find(N>5);
        %med1(jj)=plot(ax_fig(bb),BinCenter(I),Median(I),'-','linewidth',1.5,'color',colors(jj,:));
        %med2(jj)=plot(ax_fig(bb),BinCenter(I),Median(I),'ow','linewidth',0.5,'markersize',5,'markerfacecolor',colors(jj,:));
        plot(ax_fig(bb),[BinCenter(9) BinCenter(9)],[-10 170],'--k','linewidth',1.5);
        
        
        h = findobj(ax_fig(bb),'Tag','Lower Whisker');
        for j=1:length(h)-hn(bb)
            set(h(j),'Color',colors(jj,:));
        end
        h = findobj(ax_fig(bb),'Tag','Lower Adjacent Value');
        for j=1:length(h)-hn(bb)
            set(h(j),'Color',colors(jj,:));
        end    
        h = findobj(ax_fig(bb),'Tag','Upper Whisker');
        for j=1:length(h)-hn(bb)
            set(h(j),'Color',colors(jj,:));
        end
        h = findobj(ax_fig(bb),'Tag','Upper Adjacent Value');
        for j=1:length(h)-hn(bb)
            set(h(j),'Color',colors(jj,:));
        end  

        hn(bb) = numel(h);  
 
    end
    
    %uistack(med1(1),'top');
    %uistack(med2(1),'top');

end

for bb = 1:numel(basins)

    grid(ax_fig(bb),"on");
    box(ax_fig(bb),"on");
    ylim(ax_fig(bb),[0 170]);
    %xlim(ax_fig(bb),[])
    title(ax_fig(bb),glaciernames(bb));

    if bb==1
        xticklabels(ax_fig(bb),"");
        legend(ax_fig(bb),hlegend(:),["year 0","year 200"],"Location","northwest");
    elseif bb==2
        xticklabels(ax_fig(bb),"");
        yticklabels(ax_fig(bb),"");
    elseif bb==4
        yticklabels(ax_fig(bb),"");
    end

end



ylabel(tlo_fig,"Average basal met rate [m/a]","fontsize",14);
xlabel(tlo_fig,"Ice-shelf draft [m]","fontsize",14);

pos = get(H,"Position");
set(H,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);
fname = "../Figures/S3_Melt_vs_draft_Boxplot"+runID;   
print(H,fname,"-dpng","-r400");
        
