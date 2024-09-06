%% Perform Differential Analysis on RNA-Seq Data
% 

% Copyright 2020 The MathWorks, Inc.

%% 
% Use RNA-seq count data that consists of two biological
% replicates of the control (untreated) samples and two biological
% replicates of the knock-down (treated) samples [6]. Load the table with
% read counts for genes.
load("pasilla_count_noMM.mat","geneCountTable")
%% 
% Display the first few rows of the table.
head(geneCountTable,10)

%% GENERATE TABLE FOR GENE COUNTS
% miRNA = miRNAsexpressedallsamples13JUL.miRNA; 
% DS3_13JUL = miRNAsexpressedallsamples13JUL.read_count; 
% %%
% DS3_16JUN = zeros(length(miRNA),1);
% temp = miRNAsexpressedallsamples16JUN.miRNA;
% for i = 1:length(miRNA)
%     A = miRNAsexpressedallsamples16JUN.read_count(temp == miRNA(i));
%     DS3_16JUN(i) = A(1);
%     temp
% end
%%
%shave no-reads
miRNAsexpressed_13JUL(miRNAsexpressed_13JUL.read_count == 0,:) = [];
miRNAsexpressed_13JUL = sortrows(miRNAsexpressed_13JUL,"read_count");
miRNAsexpressed_16JUN(miRNAsexpressed_16JUN.read_count == 0,:) = [];
miRNAsexpressed_16JUN = sortrows(miRNAsexpressed_16JUN,"read_count");
miRNAsexpressed_20JAN(miRNAsexpressed_20JAN.read_count == 0,:) = [];
miRNAsexpressed_20JAN = sortrows(miRNAsexpressed_20JAN,"read_count");
miRNAsexpressed_30SEP(miRNAsexpressed_30SEP.read_count == 0,:) = [];
miRNAsexpressed_30SEP = sortrows(miRNAsexpressed_30SEP,"read_count");
miRNAsexpressed_DS3old(miRNAsexpressed_DS3old.read_count == 0,:) = [];
miRNAsexpressed_DS3old = sortrows(miRNAsexpressed_DS3old,"read_count");
miRNAsexpressed_hNP(miRNAsexpressed_hNP.read_count == 0,:) = [];
miRNAsexpressed_hNP = sortrows(miRNAsexpressed_hNP,"read_count");

%%
 agg = [miRNAsexpressed_13JUL.miRNA;miRNAsexpressed_16JUN.miRNA;miRNAsexpressed_20JAN.miRNA;...
     miRNAsexpressed_30SEP.miRNA;miRNAsexpressed_DS3old.miRNA;miRNAsexpressed_hNP.miRNA];
miRNA = strings(1000,1);
for i = 1:1000
    try
        [C,~,ic] = unique(agg);
        miRNA(i) = C{mode(ic)};
        agg(agg == miRNA(i)) = [];
    catch
        break
    end
end
miRNA(miRNA == "") = []; 
%% read count table
names = {"miRNAsexpressed_13JUL","miRNAsexpressed_16JUN","miRNAsexpressed_20JAN",...
    "miRNAsexpressed_30SEP","miRNAsexpressed_DS3old","miRNAsexpressed_hNP"};
A = zeros(length(miRNA),6);
%%
temp = miRNAsexpressed_hNP.miRNA;
[C,ia,ic] = unique(temp);
miRNAsexpressed_hNP = miRNAsexpressed_hNP(ia,:);
for i = 1:length(miRNA)
    readVal = miRNAsexpressed_hNP.read_count(miRNAsexpressed_hNP.miRNA == miRNA(i));
    if isempty(readVal)
        A(i,6) = 0; 
    else
        A(i,6) = readVal; 
    end
end
%%
T = table(miRNA, A(:,1), A(:,2),A(:,3),A(:,4),A(:,5));
T.Properties.VariableNames = {'ID','13JUL','16JUN','20JAN','30SEP','DS3old'};
%%
%add missing values to DS3 when compared to hNP
miRNA = setdiff(miRNAsexpressed_hNP.miRNA, miRNAsexpressed_13JUL.miRNA);
read_count = zeros(length(miRNA),1);
B = table(miRNA,read_count);
miRNAsexpressed_13JUL = [miRNAsexpressed_13JUL ; B]; 
miRNA = setdiff(miRNAsexpressed_hNP.miRNA, miRNAsexpressed_16JUN.miRNA);
read_count = zeros(length(miRNA),1);
B = table(miRNA,read_count);
miRNAsexpressed_16JUN = [miRNAsexpressed_16JUN ; B];
miRNA = setdiff(miRNAsexpressed_hNP.miRNA, miRNAsexpressed_20JAN.miRNA);
read_count = zeros(length(miRNA),1);
B = table(miRNA,read_count);
miRNAsexpressed_20JAN = [miRNAsexpressed_20JAN ; B]; 
miRNA = setdiff(miRNAsexpressed_hNP.miRNA, miRNAsexpressed_30SEP.miRNA);
read_count = zeros(length(miRNA),1);
B = table(miRNA,read_count);
miRNAsexpressed_30SEP = [miRNAsexpressed_30SEP ; B]; 
miRNA = setdiff(miRNAsexpressed_hNP.miRNA, miRNAsexpressed_DS3old.miRNA);
read_count = zeros(length(miRNA),1);
B = table(miRNA,read_count);
miRNAsexpressed_DS3old = [miRNAsexpressed_DS3old ; B]; 
%%
readcount_13JUL = miRNAsexpressedallsamples13JUL.read_count; 
readcount_16JUN = miRNAsexpressedallsamples16JUN.read_count; 
readcount_20JAN = miRNAsexpressedallsamples20JAN.read_count; 
readcount_30SEP = miRNAsexpressedallsamples30SEP.read_count; 
readcount_DS3old = miRNAsexpressedallsamplesDS3old.read_count; 
readcount_hNP = miRNAsexpressedallsampleshNP.read_count; 
%%

%% 
% Perform the differential analysis of the control and treated samples
% using the read count data for genes. Specify both replicates for each
% condition. The input |geneCountTable| has an |ID| column. Optionally, you
% can append this column to the output table by using |IDColumns| .
diffTable = rnaseqde(T,["13JUL","16JUN"],...
                     ["20JAN","30SEP","DS3old"],IDColumns="ID");
head(diffTable,5)
%%
% Look at the difference in gene expression between two conditions by
% displaying the log2 fold change for each gene.
figure
scatter(log2(mean([diffTable.Mean1,diffTable.Mean2], 2)),diffTable.Log2FoldChange,3,diffTable.AdjustedPValue,'o')
colormap(flipud(cool(256)))
colorbar;
ylabel("log2(Fold change)")
xlabel("log2(Mean of normalized counts)")
title("Fold change by FDR")
%%
% You can also annotate the values in the plot with the corresponding gene names,
% interactively select genes, and export gene lists to the workspace.
warnSettings = warning('off','bioinfo:mairplot:ZeroValues');
mairplot(diffTable.Mean2,diffTable.Mean1,Labels=T.ID,Type="MA");
set(get(gca,"Xlabel"),"String","mean of normalized counts")
set(get(gca,"Ylabel"),"String","log2(fold change)")
%%
%
warning(warnSettings);