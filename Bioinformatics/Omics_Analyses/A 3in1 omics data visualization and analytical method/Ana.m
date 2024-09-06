% 09/02/2022 the .m file 'Ana' version 1.0 is created by Hefei Zhao, PhD. This code was designed for analyzing olive pomace polyphenols concentrations from HPLC-DAD as well as many other omics data.

% Supplementary materials. S. Code 1
% A coding basis and three-in-one integrated data visualization method ‘Ana’ for the rapid analysis of multidimensional omics dataset
% Hefei Zhao1, Selina C. Wang1,*
% 1Department of Food Science and Technology, University of California, Davis, One Shields Ave, Davis, CA 95616, USA
% *Corresponding author: Selina C. Wang, email: scwang@ucdavis.edu
% First author: Hefei Zhao, email: hzhao@huskers.unl.edu; hefzhao@ucdavis.edu
% For full user guidance, please see the first reference, and it is a free/open access (OA) article.
%
% References
% Zhao, H., & Wang, S. C. (2022). A Coding Basis and Three-in-One Integrated Data Visualization Method &lsquo;Ana&rsquo; for the Rapid Analysis of Multidimensional Omics Dataset. In Life (Vol. 12, Issue 11, p. 1864). https://doi.org/10.3390/life12111864
%﻿ Zhao, H., Avena-Bustillos, R. J., & Wang, S. C. (2022). Extraction, Purification and In Vitro Antioxidant Activity Evaluation of Phenolic Compounds in California Olive Pomace. In Foods (Vol. 11, Issue 2). https://doi.org/10.3390/foods11020174
% ﻿matlab - How I obtain bars with function bar3 and different widths for each bar? - Stack Overflow. (n.d.). Retrieved September 3, 2022, from https://stackoverflow.com/questions/24269516/how-i-obtain-bars-with-function-bar3-and-different-widths-for-each-bar

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The percentage icon '%' means the contents after the % are all text notes instead of executable code.
clear, clc% clean RAM and command window
close all% close all figures, but Clustergram 1 must be closed manully.
tic % start timing
%%%%%%%%%%%%%%%%%%
% Initial input area for normal users or beginners

%fn= 'ab126omics_forAna';% input the excel file name in the '' area, olivephenolics can be replaced by user's excel file name, the excel file must be in the same folder as the this .m file
%fn = 'ab126omics_forAna_normC';
fn = 'ab126omics_forAna_threshold_newMS';
unt= 'Score';% input unit of data in the '' area, mg/g can be replaced by the user's unit, such as %, g/mL, mg/mL etc.
fs= 10;% input font size, 18-22 are recommended, must be >=7
cl= parula;% color style: jet(256) is rainbow; cool is blue to pink;  parula is blue to yellow; redbluecmap is blue to red; [] is transparent; user can replaced jet(256) by cool, [] or parula or bredbluecmap, etc.
pcamz= 10;% PCA marker size
pcalable= 0;% 0 will lable PCA vectors by variable/ compound numbers; 1 will lable PCA vectors by variable/ compound full names.
mk= '.'; % set PCA data marker in '' area; . is dot; p is star; s is square; * is snow flasker; o is o; 'd' is rhombus.

% Note: Method for print or output figures:

% Figure 1, 3D heatmap
% print(fig3Dbar,'olive-adjudtsize.png','-dpng','-r150');% -r150 defines 150dpi,-r300 will provide 300 dpi, -r100 will provide 100 dpi, etc. See S. Fig. 3

% Clustergram 1, heatmap cluster, can only be output to a pdf for high-resolution as descripted in S. Fig. 4.

% Figure 2-7, print PCA charts automatically at 150 dpi

%%%%%%%%%%%%%%%%%%% 
% It may result in errors if modify code below this line if not skillful on MATLAB.
% But it would be greatly helpful for active learners to modify the code to improve coding skills.

fname= strcat(fn);% excel file name
info= readtable(fname,'PreserveVariableNames',true, 'ReadRowNames',true);% read data again for column names
rm= string(info.Properties.RowNames);% read row names of samples
cm= string(info.Properties.VariableNames);% read column/variable names of phenolic compounds
sz= size(info);% read data area size
ns= sz(1);% read sample number
nb= sz(2);% read compound number
num= readmatrix(fname); % read data only
num(:,1) = []; % delete first column NAN
ave= num; % define ave values
%ave = normalize(ave); 
%

% plot 3D bar chart with size adjustments
%
fig3Dbar=figure; % name figure as fig3Dbar
h=bar3(ave,'detached');%'detached', 'grouped', or 'stacked'
%
values=ave; % values equal to ave.
m= 1.1*max(values(:))*2; % normalize constant for bar width, 1.1 ensure a small distance between bars
shading interp % set color shading properties
for i = 1:length(h)% the following code set bottom size in accordance with the z-axis values
    % Get the ZData matrix of the current group
    xdata= get(h(i),'Xdata');
    ydata= get(h(i),'Ydata');
    zdata= get(h(i),'Zdata');
    set(h(i),'Cdata',zdata)
    for k = 1:size(xdata,1)/6 
        datax = xdata((k-1)*6+1+(0:5),:);
        datay = ydata((k-1)*6+1+(0:5),:);
        dirx=((datax-round(datax))>0)-((datax-round(datax))<0);
        diry=((datay-round(datay))>0)-((datay-round(datay))<0);
        xdata((k-1)*6+1+(0:5),:) = round(xdata((k-1)*6+1+(0:5),:),0)+dirx*values(ceil(((k-1)*6+1)/6),i)/m;
        ydata((k-1)*6+1+(0:5),:) = round(ydata((k-1)*6+1+(0:5),:),0)+diry*values(ceil(((k-1)*6+1)/6),i)/m;
    end
    set(h(i),'XData',xdata);
    set(h(i),'YData',ydata);
end
%set(h,'EdgeColor','k') % set edge color as k black
view(-25, 83); % default view angle
colormap (cl) % set colormap as cl
colorbar % show color bar
%
xticks(1:nb)% add x ticks 
yticks(1:ns)% add y ticks 
set(gca,'XTickLabel',cm)% label x axis by column/variable names of phenolic compound
set(gca,'YTickLabel',rm)% label y axis by row names of samples
set(h,'FaceAlpha',.8) % set transparency of bars to 0.5
zlabel(unt,'FontSize',fs) % set font size of z-axis
ax = gca;
ax.FontSize = fs; % set font size of color code bar
cb=colorbar;
colormap(cl); % define color as descripted by cl
cb.Label.String = unt; % set unit of color code bar as descripted by unt
%print(fig3Dbar,'olive-adjudtsize.png','-dpng','-r150');% -r150 defines 150dpi,-r300 will provide 300 dpi, -r100 will provide 100 dpi, etc.

% Cluster analysis
cfac=clustergram(ave, 'RowLabels', rm,'ColumnLabels',cm,'Colormap',colormap(cl),'Standardize','Row'); % the code standardize data on each row
% For more information see: https://www.mathworks.com/help/bioinfo/ref/clustergram.html
set(cfac,'Linkage','complete','Dendrogram',5)
set(cfac,'Annotate','off') % turn of annotate
set(cfac,'Linkage','Average') % set linkage method by Average
set(cfac,'RowPDist','Euclidean') % row distance method Euclidean
set(cfac,'ColumnPDist','Euclidean') % column distance method Euclidean
%
% PCA
sdz = zscore(ave,[],2); % standardized data along data rows
[coefs,score,latent,tsquared,explainedvariance] = pca(sdz); % run PCA
%
if pcalable== 0 % 0 will lable PCA vectors by variable/ compound numbers; 1 will lable PCA vectors by variable/ compound full names.
    lbls= string(1:length(cm));
else
    lbls= cm;
end
%
figPCA12biplot= figure; % open a new figure
pa12= biplot(coefs(:,1:2),'Scores',score(:,1:2),'VarLabels',lbls,'Marker',mk,'MarkerSize',pcamz); % plot PCA biplot of PC1 vs. PC2
grid off % turn off grid
box off % ture off box
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel(['PC1',' ', num2str(explainedvariance(1)),'%'],'FontSize',fs)% label x-axis
ylabel(['PC2',' ', num2str(explainedvariance(2)),'%'],'FontSize',fs)% label y-axis
print(figPCA12biplot,'PCA12-biplot.png','-dpng','-r150');% print figure at 150 dpi
%
clr = hsv(ns);
figPCA12score= figure;
gscatter(score(:,1),score(:,2),rm,clr,mk) % plot PCA scoreplot of PC1 vs. PC2
legend('Location','northeastoutside')
grid off % turn off grid
box off % ture off box
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel(['PC1',' ', num2str(explainedvariance(1)),'%'],'FontSize',fs)% label x-axis
ylabel(['PC2',' ', num2str(explainedvariance(2)),'%'],'FontSize',fs)% label y-axis
print(figPCA12score,'PCA12-score.png','-dpng','-r150');% print figure at 150 dpi
%
figPCA23biplot= figure; % open a new figure
pa23= biplot(coefs(:,2:3),'Scores',score(:,2:3),'VarLabels',lbls,'Marker',mk,'MarkerSize',pcamz); % plot PCA biplot of PC2 vs. PC3
grid off % turn off grid
box off % ture off box
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel(['PC2',' ', num2str(explainedvariance(2)),'%'],'FontSize',fs)
ylabel(['PC3',' ', num2str(explainedvariance(3)),'%'],'FontSize',fs)
print(figPCA23biplot,'PCA23-biplot.png','-dpng','-r150');
%
figPCA12score= figure;
gscatter(score(:,2),score(:,3),rm,clr,mk)% plot PCA scoreplot of PC2 vs. PC3
legend('Location','northeastoutside')
grid off % turn off grid
box off % ture off box
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel(['PC2',' ', num2str(explainedvariance(2)),'%'],'FontSize',fs)% label x-axis
ylabel(['PC3',' ', num2str(explainedvariance(3)),'%'],'FontSize',fs)% label y-axis
print(figPCA12score,'PCA23-score.png','-dpng','-r150');% print figure at 150 dpi
%
figPCA123biplot= figure; % open a new figure
pa123= biplot(coefs(:,1:3),'Scores',score(:,1:3),'VarLabels',lbls,'Marker',mk,'MarkerSize',pcamz); % plot PCA biplot of PC1 vs. PC2 vs. PC3
grid off % turn off the grid
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel(['PC1',' ', num2str(explainedvariance(1)),'%'],'FontSize',fs)
ylabel(['PC2',' ', num2str(explainedvariance(2)),'%'],'FontSize',fs)
zlabel(['PC3',' ', num2str(explainedvariance(3)),'%'],'FontSize',fs)% label z-axis
box on % ture on box
print(figPCA123biplot,'PCA123-biplot.png','-dpng','-r150');
%
ev= figure;
pareto(explainedvariance,0.99) % draw Explained Variance
ax = gca;
ax.FontSize = fs-6; % set font size of tick, 6 less than label font size
xlabel('Principal Component','FontSize',fs)
ylabel('Explained Variance(%)','FontSize',fs)
box off % ture off box
print(ev,'PCA-Explained-Variance.png','-dpng','-r150');
toc% stop timing