%% Ridgeline Testing
%joyPlot(data',x,4,'FaceColor',mean(data,"all"))
%colorbar
%%
%table = readtable("Ridgeline Test.xlsx");
%table = readtable("Bin Data DP15AUG2023.xlsx");
table = readtable("Bin Data DP07MAY2024.xlsx");
table = readtable("DP02JUL2024 FT Bin Data.xlsx");
names = unique(table.Sample);
data = nan(1000,numel(names));
areas = nan(numel(names),1);
for i = 1 : numel(names)
    if contains(names{i},'FT') %&& contains(names{i},'0M')
        v = table(table.Sample == string(names{i}),'ConcentrationAverage');
        data(:,i) = v{:,:};

        bin = table(1:1000,"BinCentre_nm_");
        binx = bin{:,:};
        areas(i,1) = polyarea(binx,data(:,i));
    end
end
%%
areas = areas(all(~isnan(areas),2),:);
tickNames = names(contains(names,'FT')); %& contains(names,'0M'));
data = data(:,all(~isnan(data)));
data = data(51:500,:);
%%
%x = linspace(0,93,size(data,2));
bin = table(51:500,"BinCentre_nm_");
binx = bin{:,:};
x = linspace(50.5,499.5,450);
%areas = polyarea(repmat(binx,1,size(data,2)),data);
%[hs,hf] = joyPlot(flip(data,2),x,0.01,'FaceColor',flip(areas,1),'LineColor','k');

%color plots by sample ID
%[hs,hf] = joyPlot(data,x,0.01,'FaceColor',normalize(mean(data,1),'range'),'LineColor','k');
[hs,hf] = joyPlot(data,x,0.01,'FaceColor',areas,'LineColor','k');
%[hs,hf] = joyPlot(data,x,0.01,'FaceColor','position','LineColor','k');
%legend(tickNames)
%tickNames = flip(tickNames,1);
for i = 1:numel(hs)
    hs(i).DisplayName = tickNames{i};
end
colormap('parula')
legend(hs)
cbar = colorbar;
%set(gcf,'position',[500,100,560,680],'InvertHardcopy','off')
set(gcf,'Color', 'none')
set(gca,'Color','none')
set(gca,'XGrid','on')
set(gca,'YTick',[])
set(hs,'FaceAlpha',0.3)
%set(gca,'Visible','off','box','off', 'XTick',binx,'YTick',[])