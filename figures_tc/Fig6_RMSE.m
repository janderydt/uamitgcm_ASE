
function Fig5_RMSE

RMSE=[0.10,.340,.10,.060,.230,.240,.070,.08;...
0.210,.150,.10,.080,.190,.250,.220,.22;...
0.060,.080,.020,.030,.270,.170,.070,.06;...
0.40,.280,.080,.10,.440,.680,.970,.54]';


RMSE_hi = RMSE(1:2:7,:);
RMSE_av = RMSE(2:2:8,:);

Hfig = fig('units','inches','height',30.84*12/72.27,'width',70.1*12/72.27,'fontsize',12,'font','Helvetica');

tlo_fig = tiledlayout(1,2,"TileSpacing","compact");
for i = 1:2
    ax_fig(i) = nexttile(tlo_fig,i); hold on;
end

b=bar(ax_fig(1),RMSE_hi); title(ax_fig(1),"$\mbox{\textit{hi\_melt}}$","interpreter","latex")
bar(ax_fig(2),RMSE_av); title(ax_fig(2),"$\mbox{\textit{av\_melt}}$","interpreter","latex")

for i=1:2
    grid(ax_fig(i),"on");
    box(ax_fig(i),"on");
    xticks(ax_fig(i),[1:4])
    xticklabels(ax_fig(i),{'Pine Island','Thwaites','Crosson','Dotson'})
end
yticklabels(ax_fig(2),'')
ylabel(ax_fig(1),"Root Mean Square Error");
lg=legend(b,["$\;{T}_{\star {\rm IF}}\rightarrow\,{T}_{\star {\rm IF}}(t=0)\quad$","$\;\mu\rightarrow\,\mu(t=0)\quad$","$\;\epsilon_T\rightarrow\,\epsilon_T(t=0)\quad$","$\;\epsilon_U\rightarrow\,\epsilon_U(t=0)$"],"interpreter","latex","Orientation","horizontal","FontSize",14);
lg.Layout.Tile="north";

pos = get(Hfig,"Position");
set(Hfig,"PaperPositionMode","Auto","PaperUnits","Inches","PaperSize",[pos(3),pos(4)]);

fname = "../Figures/Fig5_RMSE";
print(Hfig,fname,"-dpng","-r400");

