%% Prepare Data from 09JAN2023
path_09JAN2023 = "C:\Users\StevenSummey\Documents\Characterization\NTA\2023\JANUARY\09JAN2023\Full Plate-ExperimentSummary_09JAN2023.xlsx";
sampleNames_09JAN2023 = table2array(readtable(path_09JAN2023,"Range","B12:AO12","ReadVariableNames",false));
sampleNames_09JAN2023 = sampleNames_09JAN2023(:);
sampleTotalConcentration_09JAN2023 = table2array(readtable(path_09JAN2023, "Range","B17:AO17"));
sampleTotalConcentration_09JAN2023 = sampleTotalConcentration_09JAN2023(:); 
expectedTotalConcentration_09JAN2023 = ones(numel(sampleTotalConcentration_09JAN2023),1).*(1.33E12); %expected concentration from manual 
                                                                                                      %measurement of DS3-12DEC2022B-1
dilutionFactor = ones(numel(sampleTotalConcentration_09JAN2023),1);                                                                                                      
for i = 1:numel(sampleTotalConcentration_09JAN2023)
    dilutionString = strsplit(sampleNames_09JAN2023{i},' ');
    dilutionString = regexp(dilutionString{2},'\d+','match');
    dilutionFactor(i,1) = cellfun(@str2num,dilutionString);
end
sampleLoadedConcentration_09JAN2023 = sampleTotalConcentration_09JAN2023./dilutionFactor;
expectedLoadedConcentration_09JAN2023 = expectedTotalConcentration_09JAN2023./dilutionFactor; 
invDilutionFactor_09JAN2023 = 1./dilutionFactor;

% Parity Plot
%CREATEFIT(EXPECTEDLOADEDCONCENTRATION_09JAN2023,SAMPLELOADEDCONCENTRATION_09JAN2023)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: expectedLoadedConcentration_09JAN2023
%      Y Output: sampleLoadedConcentration_09JAN2023
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Jan-2023 14:26:45


%% Fit: 'Parity'.
[xData, yData] = prepareCurveData( expectedLoadedConcentration_09JAN2023, sampleLoadedConcentration_09JAN2023 );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Create a figure for the plots.
figure( 'Name', 'Parity Plot' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData);
legend( h, 'Sample Loaded Conc, (p/mL) vs. Expected Loaded Conc, (p/mL)', 'Parity', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Expected Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'Sample Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'Parity - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Expected Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'Sample Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on
%--------------------------------------------------------------------------------------------------------------%
%% Manual Data from two (2) operators, include expected
%AB126 Dilution Factors
%1000, 2500, 5000, 10000, 25000
dfAB126_manual = manualOperator1_Dilutions;
invDFAB126 = 1./dfAB126_manual; 
%expected concentrations
expConcAB126 = 1.49E12./dfAB126_manual; %ds3-08SEP2022B-1
%loaded particle concentration
%ALP
loadedConcALP = manualOperator2_loadedParticles;
%SS
%loadedConcSS = manualOperator1_loadedParticles;

%from plate
loadedConcPlate = manualOperatorFromPlate_loadedParticles;

combinedData = [loadedConcALP, loadedConcSS];
dimData = size(combinedData);
numColumns = dimData(2);
%%
xData = invDFAB126;
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

plot(xData, expConcAB126,'*','Color','k');
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
legend(legendCell, 'Location', 'northwest');


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


%-------------------------------------------------------------------------------------------------%
%% Manual Run Sampled from Plate compared to AutoSampler Data, include expected
%PARITY PLOT
expectedTotalConcentration_05JAN2023 = ones(numel(autoSampler05JAN2023_totalConcentration),1).*(1.33E12);
expectedLoadedConcentration_05JAN2023 = expectedTotalConcentration_05JAN2023./autoSampler05JAN2023_Dilutions; 
invDilutionFactor_05JAN2023 = 1./autoSampler05JAN2023_Dilutions;
% Fit: 'Parity'.
[xData, yData] = prepareCurveData( manualOperatorFromPlate_loadedParticles, autoSampler05JAN2023_loadedParticles );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Create a figure for the plots.
figure( 'Name', 'Parity Plot - Plate Samples' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData);
legend( h, 'AutoSampler Loaded Conc, (p/mL) vs. Manual Loaded Conc, (p/mL)', 'Parity', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Manual Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'AutoSampler Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'Parity - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Manual Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'AutoSampler Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on
%--------------------------------------------------------------------------------------------------------------%




%% New AutoSampler Run w/ EDTA compared to expected particle concentration and other manual runs
%PARITY PLOT
path_19JAN2023 = "C:\Users\StevenSummey\Documents\Characterization\NTA\2023\JANUARY\19JAN2023\Full Plate-ExperimentSummary_19JAN2023.xlsx";
% Prepare Data from 09JAN2023
%path_09JAN2023 = "C:\Users\StevenSummey\Documents\Characterization\NTA\2023\JANUARY\09JAN2023\Full Plate-ExperimentSummary_09JAN2023.xlsx";
sampleNames_19JAN2023 = table2array(readtable(path_19JAN2023,"Range","B14:AT14","ReadVariableNames",false));
sampleNames_19JAN2023 = sampleNames_19JAN2023(:);
sampleLoadedConcentration_19JAN2023 = table2array(readtable(path_19JAN2023, "Range","B16:AT16"));
sampleLoadedConcentration_19JAN2023 = sampleLoadedConcentration_19JAN2023(:); 
expectedTotalConcentration_19JAN2023 = ones(numel(sampleLoadedConcentration_19JAN2023),1).*(1.33E12); %expected concentration from manual 
                                                                                                      %measurement of DS3-12DEC2022B-1
dilutionFactor = ones(numel(sampleLoadedConcentration_19JAN2023),1);                                                                                                      
for i = 1:numel(sampleLoadedConcentration_19JAN2023)
    dilutionString = strsplit(sampleNames_19JAN2023{i},' ');
    dilutionString = regexp(dilutionString{2},'\d+','match');
    dilutionFactor(i,1) = cellfun(@str2num,dilutionString);
end
sampleTotalConcentration_19JAN2023 = sampleLoadedConcentration_19JAN2023.*dilutionFactor;
expectedLoadedConcentration_19JAN2023 = expectedTotalConcentration_19JAN2023./dilutionFactor; 
invDilutionFactor_19JAN2023 = 1./dilutionFactor;
% Parity Plot
%CREATEFIT(EXPECTEDLOADEDCONCENTRATION_09JAN2023,SAMPLELOADEDCONCENTRATION_09JAN2023)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: expectedLoadedConcentration_09JAN2023
%      Y Output: sampleLoadedConcentration_09JAN2023
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Jan-2023 14:26:45


% Fit: 'Parity'.
[xData, yData] = prepareCurveData( expectedLoadedConcentration_19JAN2023, sampleLoadedConcentration_19JAN2023 );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Create a figure for the plots.
figure( 'Name', 'Parity Plot' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult, xData, yData);
legend( h, 'Sample Loaded Conc, (p/mL) vs. Expected Loaded Conc, (p/mL)', 'Parity', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Expected Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'Sample Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult, xData, yData, 'residuals' );
legend( h, 'untitled fit 1 - residuals', 'Zero Line', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Expected Loaded Conc, (p/mL)', 'Interpreter', 'none' );
ylabel( 'Sample Loaded Conc, (p/mL)', 'Interpreter', 'none' );
grid on