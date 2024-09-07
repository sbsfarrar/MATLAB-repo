%100nm Dilution Factors
%1E7, 5E7, 7.5E7, 1.5E78, 3E8, 6E8, 1E9
df100nm = [1818914; 363783; 242522; 121261; 60630; 30315; 18189];
invDF100nm = 1./df100nm; 
%expected concentrations
expConc100nm = 1.82E13./df100nm; 
%loaded particle concentration
%ALP
loadedConcALP = [1.67E7; 5.62E7; 1.34E8; 2.41E8; 4.67E8; 8.97E8; 1.30E9];
%SS
loadedConcSS = [1.14E7; 5.84E7; 8.51E7; 1.68E8; 3.26E8; 5.25E8; 9.41E8];

combinedData = [loadedConcALP, loadedConcSS];
dimData = size(combinedData);
numColumns = dimData(2);
%%
xData = invDF100nm;
legendCell = cell(numColumns*2,1);
eqnTextCell = cell(numColumns,1);
fitresultStore = cell(numColumns,1);
for i = 1:numColumns
    yData = combinedData(:,i);
    [fitresult, gof] = createFit_Concentration(xData,yData);
    fitresultStore{i} = fitresult;
    p = plot(fitresult);
    hold on
    if i == 1
        set(p,'Color','r');
        plot(xData,yData,'o','Color','r');
        legendCell{1} = ['Analyst 1, R-squared = ' num2str(gof.rsquare)];
        legendCell{2} = 'Analyst 1 Data';
        fit1 = fitresult(1);
        fit2 = fitresult(2);
        eqnTextCell{1} = sprintf('y = %f*x + %f',string(num2str(fit1,'%.e ')),string(num2str(fit2,'%.e ')));
    elseif i == 2
        set(p,'Color','g');
        plot(xData,yData,'o','Color','g');
        legendCell{3} = ['Analyst 2, R-squared = ' num2str(gof.rsquare)];
        legendCell{4} = 'Analyst 2 Data';
        fit1 = fitresult(1);
        fit2 = fitresult(2);
        eqnTextCell{2} = sprintf('y = %f*x + %f',string(num2str(fit1,'%.e ')),string(num2str(fit2,'%.e ')));
    else
        set(p,'Color','b');
        plot(xData,yData,'o','Color','b');
        legendCell{5} = ['Analyst 3, R-squared = ' num2str(gof.rsquare)];
        legendCell{6} = 'Analyst 3 Data';
        fit1 = fitresult(1);
        fit2 = fitresult(2);
        eqnTextCell{3} = sprintf('y = %f*x + %f',string(num2str(fit1,'%.e ')),string(num2str(fit2,'%.e ')));
    end
end

plot(xData, expConc100nm,'*','Color','k');
legendCell{end+1} = 'Expected Particle Concentration';

% %Add equation text
% xl = xlim;
% yl = ylim;
% xt = 0.05*(xl(2)-xl(1)) + xl(1);
% yt = 0.90*(yl(2)-yl(1)) + yl(1);
% t = text(xt,yt,eqnTextCell,"FontSize",11);
% for i = 1:numColumns
%     if i == 1
%         t(i).Color = 'r';
%     elseif i == 2
%         t(i).Color = 'g';
%     else
%         t(i).Color = 'b';
%     end
% end
% set(gca);

title('Size Standard Concentration Standard Curve');
xlabel('1 / Dilution Factor');
ylabel('Measured Loaded Particles (particles/mL)');
legend(legendCell, 'Location', 'bestoutside');


% for i = 1:numColumns
%     yData = combinedData(:,i);
%     if i == 1
%         p1 = plot(xData,yData,'o','Color','r');
%         plotArray = p1;
%     elseif i == 2
%         p2 = plot(xData,yData,'o','Color','g');
%         plotArray = [plotArray p2];
%     else
%         p3 = plot(xData,yData,'o','Color','b');
%         plotArray = [plotArray p3];
%     end
% end
% a = axes('Position',get(gca,'Position'),'Visible','off');
% legend(a,plotArray,eqnTextCell,'Location','bestoutside');
hold off
%ADD CREATE FIT FUNCTION
%%
h1 = plot(fittedmodel_ALP); % h1 contains the lines for model 1
set(h1,'Color','r'); % make data line for model 1 green
%set(h1(2:end),'Color','c'); % make fit and confidence lines for model 1 cyan
hold on
h2 = plot(fittedmodel_SS); % h2 contains the lines for model 2
set(h2,'Color','g'); % make data line for model 2 orange
%set(h2(2:end),'Color','m'); % make fit and confidence lines for model 2 magenta
legend({['Analyst 1, R-square = ' num2str(goodness_ALP.rsquare)],...
    ['Analyst 2, R-square = ' num2str(goodness_SS.rsquare)]},'Location','northwest');
