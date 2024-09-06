function plotCounts(fcCountsTable,cufflinksGenesFPKMFile)
    genesFPKMTable = readtable(cufflinksGenesFPKMFile,FileType="text");
    % Plot counts of genes identified by Cufflinks.
    figure
    geneNames = categorical(fcCountsTable.ID,fcCountsTable.ID);
    stem(geneNames, log2(fcCountsTable.Aligned_sorted))
    xlabel("Cufflinks-identified genes")
    ylabel("log2 counts")
    
    % Plot counts along their respective genomic positions.
    geneStart = str2double(extractBetween(genesFPKMTable.locus,":","-"));
    figure
    stem(geneStart,log2(fcCountsTable.Aligned_sorted))
    xlabel("Reference Genome Position")
    ylabel("log2 counts")
end