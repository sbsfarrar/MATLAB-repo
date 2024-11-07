%% LOAD IMAGES AND SET UP STRUCTURE
%set folder paths
%\\10.10.110.10\NTAData\Quality Control\2024-10-14-SS
dirPathQC = "\\10.10.110.10\NTAData\Quality Control"; 
dirPathPD = "\\10.10.110.10\NTAData\Process Development";
dirPathArch = "\\10.10.110.10\NTAData\Quality Control\#ARCHIVE";
imagesQC = dir(dirPathQC + "\*\*background*.bmp");
imagesPD = dir(dirPathPD + "\*\*background*.bmp");
imagesArch = dir(dirPathArch + "\*\*background*.bmp");
images = [imagesQC; imagesPD; imagesArch];

trackingStruct = struct();
%parse operater, date, and sample name
for i=1:length(images)
    disp(i)
    try
        folderParse = strsplit(images(i).folder,"\");
        runString = folderParse{end};
        operator = string(regexp(runString,'[A-Z]+', 'match'));
        date = string(regexp(runString,'[0-9]+-[0-9]+-[0-9]+','match'));
        sample_summ = strrep(images(i).name,"-background.bmp","");
        sample_summ = strrep(sample_summ,regexp(sample_summ,"[0-9]+-[0-9]+-[0-9]+ [0-9]+-[0-9]+-[0-9]+","match"),"");
        sample_summ = strtrim(sample_summ);
        imagePath = images(i).folder + "\" + images(i).name;
        trackingStruct.Sample(i,1) = sample_summ; 
        trackingStruct.Date(i,1) = date; 
        trackingStruct.Operator(i,1) = operator; 
        trackingStruct.Path(i,1) = imagePath;
    catch
        %remove problematic file
        trackingStruct.Sample(i,1) = missing; 
        trackingStruct.Date(i,1) = missing; 
        trackingStruct.Operator(i,1) = missing; 
        trackingStruct.Path(i,1) = missing;
        disp("skipping")
        continue;
    end
end

disp("removing missing vals")
names = fieldnames(trackingStruct);
for i=1:length(names)
    trackingStruct.(names{i}) = rmmissing(trackingStruct.(names{i}));
end
%% AFTER INITIAL READ
% Filter out files already processed! (see: TFF file parse)
% read backgroud file into workspace
nta_path = "backgroundTracking.xlsx";
nta_particle_data = readtable(nta_path);

files_to_process = setdiff(trackingStruct.Path, nta_particle_data.Path);

%filter out files from "trackingStruct" to pass into next section
trackingStruct_filt = struct(); 
counter = 0; 
for i=1:length(trackingStruct.Path)
    for j=1:length(files_to_process)
        if trackingStruct.Path(i) == files_to_process(j)
            counter = counter + 1; %increment hits for indexing
            trackingStruct_filt.Sample(counter,1) = trackingStruct.Sample(i);
            trackingStruct_filt.Date(counter,1) = trackingStruct.Date(i); 
            trackingStruct_filt.Operator(counter,1) = trackingStruct.Operator(i); 
            trackingStruct_filt.Path(counter,1) = trackingStruct.Path(i);
        end
    end
end
%% READ IMAGES, SUMMARY, AND PARTICLE DATA
summaryQC = dir(dirPathQC + "\*\*_Summary.csv*");
summaryPD = dir(dirPathPD + "\*\*_Summary.csv*");
summaryArch = dir(dirPathArch + "\*\*_Summary.csv*");
summary = [summaryQC; summaryPD; summaryArch];
%particle data
particlesQC = dir(dirPathQC + "\*\*_ParticleData.csv*");
particlesPD = dir(dirPathPD + "\*\*_ParticleData.csv*");
particlesArch = dir(dirPathArch + "\*\*_ParticleData.csv*");
particles = [particlesQC; particlesPD; particlesArch];

grayAVG = zeros(length(trackingStruct_filt.Path),1);
for i=1:length(trackingStruct_filt.Path)
    disp(i)
    image = imread(trackingStruct_filt.Path(i));
    %[grayAVG(i,1), grayAVG(i,2)] = background_vals(image);
    [trackingStruct_filt.GrayVal(i,1)] = background_vals(image);
    for j=1:length(summary)
        %disp(length(summary) - 1);
        summaryPath = summary(j).folder + "\" + summary(j).name;
        sample_summ = strrep(summary(j).name,"_Summary.csv","");
        sample_summ = strtrim(sample_summ);
        if contains(trackingStruct_filt.Path(i),sample_summ)
            try
                trackingStruct_filt.Concentration(i,1) = table2array(readtable(summaryPath,"Range","B46:B46"));
                trackingStruct_filt.PPF(i,1) = table2array(readtable(summaryPath,"Range","B48:B48"));
                trackingStruct_filt.CompletedTracks(i,1) = table2array(readtable(summaryPath,"Range","B50:B50"));
                trackingStruct_filt.xDrift(i,1) = table2array(readtable(summaryPath,"Range","B51:B51"));
                trackingStruct_filt.yDrift(i,1) = table2array(readtable(summaryPath,"Range","B52:B52"));
                trackingStruct_filt.Mean(i,1) = table2array(readtable(summaryPath,"Range","B81:B81"));
                trackingStruct_filt.Mode(i,1) = table2array(readtable(summaryPath,"Range","B82:B82"));
                trackingStruct_filt.SD(i,1) = table2array(readtable(summaryPath,"Range","B83:B83"));
                trackingStruct_filt.d10(i,1) = table2array(readtable(summaryPath,"Range","B84:B84"));
                trackingStruct_filt.d50(i,1) = table2array(readtable(summaryPath,"Range","B85:B85"));
                trackingStruct_filt.d90(i,1) = table2array(readtable(summaryPath,"Range","B86:B86"));
                trackingStruct_filt.ValidTracks(i,1) = table2array(readtable(summaryPath,"Range","B87:B87"));
                
                disp("Summary Loaded")
            catch
                trackingStruct_filt.Concentration(i,1) = missing;
                trackingStruct_filt.PPF(i,1) = missing;
                trackingStruct_filt.CompletedTracks(i,1) = missing;
                trackingStruct_filt.xDrift(i,1) = missing;
                trackingStruct_filt.yDrift(i,1) = missing;
                trackingStruct_filt.Mean(i,1) = missing;
                trackingStruct_filt.Mode(i,1) = missing;
                trackingStruct_filt.SD(i,1) = missing;
                trackingStruct_filt.d10(i,1) = missing;
                trackingStruct_filt.d50(i,1) = missing;
                trackingStruct_filt.d90(i,1) = missing;
                trackingStruct_filt.ValidTracks(i,1) = missing;
                disp("skipping summary")
            end
            break;
        else
            trackingStruct_filt.Concentration(i,1) = 999999999;
            trackingStruct_filt.PPF(i,1) = 999999999;
            trackingStruct_filt.CompletedTracks(i,1) = 999999999;
            trackingStruct_filt.xDrift(i,1) = 999999999;
            trackingStruct_filt.yDrift(i,1) = 999999999;
            trackingStruct_filt.Mean(i,1) = 999999999;
            trackingStruct_filt.Mode(i,1) = 999999999;
            trackingStruct_filt.SD(i,1) = 999999999;
            trackingStruct_filt.d10(i,1) = 999999999;
            trackingStruct_filt.d50(i,1) = 999999999;
            trackingStruct_filt.d90(i,1) = 999999999;
            trackingStruct_filt.ValidTracks(i,1) = 999999999;
        end
    end
    for k=1:length(particles)
        particlePath = particles(k).folder + "\" + particles(k).name;
        sample_part = strrep(particles(k).name,"_ParticleData.csv","");
        sample_part = strtrim(sample_part);
        if contains(trackingStruct_filt.Path(i), sample_part)
            try
                particleTable = readtable(particlePath);
                particleTable_FALSE = particleTable(string(particleTable.IncludedInDistribution_) == 'False', :);
                particleTable_TRUE = particleTable(string(particleTable.IncludedInDistribution_) == 'True', :);
                trackingStruct_filt.AvgDiffCoeff_ALL(i,1) = mean(particleTable.DiffusionCoefficient_nm_2S__1);
                trackingStruct_filt.AvgDiffCoeff_FALSE(i,1) = mean(particleTable_FALSE.DiffusionCoefficient_nm_2S__1);
                trackingStruct_filt.AvgDiffCoeff_TRUE(i,1) = mean(particleTable_TRUE.DiffusionCoefficient_nm_2S__1);
                trackingStruct_filt.AvgLogIntensity_ALL(i,1) = mean(particleTable.Ln_AdjustedIntensity__AU);
                trackingStruct_filt.AvgLogIntensity_FALSE(i,1) = mean(particleTable_FALSE.Ln_AdjustedIntensity__AU);
                trackingStruct_filt.AvgLogIntensity_TRUE(i,1) = mean(particleTable_TRUE.Ln_AdjustedIntensity__AU);
                trackingStruct_filt.AvgTrackLength_ALL(i,1) = mean(particleTable.Tracklength);
                trackingStruct_filt.AvgTrackLength_FALSE(i,1) = mean(particleTable_FALSE.Tracklength);
                trackingStruct_filt.AvgTrackLength_TRUE(i,1) = mean(particleTable_TRUE.Tracklength);
                trackingStruct_filt.AvgSize_ALL(i,1) = mean(particleTable.Size_nm);
                trackingStruct_filt.AvgSize_FALSE(i,1) = mean(particleTable_FALSE.Size_nm);
                trackingStruct_filt.AvgSize_TRUE(i,1) = mean(particleTable_TRUE.Size_nm);
                
                disp("Particle Data Loaded")
            catch
                trackingStruct_filt.AvgDiffCoeff_ALL(i,1) = missing;
                trackingStruct_filt.AvgDiffCoeff_FALSE(i,1) = missing;
                trackingStruct_filt.AvgDiffCoeff_TRUE(i,1) = missing;
                trackingStruct_filt.AvgLogIntensity_ALL(i,1) = missing;
                trackingStruct_filt.AvgLogIntensity_FALSE(i,1) = missing;
                trackingStruct_filt.AvgLogIntensity_TRUE(i,1) = missing;
                trackingStruct_filt.AvgTrackLength_ALL(i,1) = missing;
                trackingStruct_filt.AvgTrackLength_FALSE(i,1) = missing;
                trackingStruct_filt.AvgTrackLength_TRUE(i,1) = missing;
                trackingStruct_filt.AvgSize_ALL(i,1) = missing;
                trackingStruct_filt.AvgSize_FALSE(i,1) = missing;
                trackingStruct_filt.AvgSize_TRUE(i,1) = missing;

                disp("skipping particles")
            end
            break;
        else
            trackingStruct_filt.AvgDiffCoeff_ALL(i,1) = 999999999;
            trackingStruct_filt.AvgDiffCoeff_FALSE(i,1) = 999999999;
            trackingStruct_filt.AvgDiffCoeff_TRUE(i,1) = 999999999;
            trackingStruct_filt.AvgLogIntensity_ALL(i,1) = 999999999;
            trackingStruct_filt.AvgLogIntensity_FALSE(i,1) = 999999999;
            trackingStruct_filt.AvgLogIntensity_TRUE(i,1) = 999999999;
            trackingStruct_filt.AvgTrackLength_ALL(i,1) = 999999999;
            trackingStruct_filt.AvgTrackLength_FALSE(i,1) = 999999999;
            trackingStruct_filt.AvgTrackLength_TRUE(i,1) = 999999999;
            trackingStruct_filt.AvgSize_ALL(i,1) = 999999999;
            trackingStruct_filt.AvgSize_FALSE(i,1) = 999999999;
            trackingStruct_filt.AvgSize_TRUE(i,1) = 999999999;
        end
    end
    disp("-----------------------------------------------------------------------")
end

%% FINAL CLEAN AND EXPORT
%fixed later, but transpose:
trackingStruct2 = trackingStruct_filt;

% trackingStruct2.Concentration = trackingStruct2.Concentration.';
% trackingStruct2.PPF = trackingStruct2.PPF.';
% trackingStruct2.CompletedTracks = trackingStruct2.CompletedTracks.';
% trackingStruct2.xDrift = trackingStruct2.xDrift.';
% trackingStruct2.yDrift = trackingStruct2.yDrift.';
% trackingStruct2.Mean = trackingStruct2.Mean.';
% trackingStruct2.Mode = trackingStruct2.Mode.';
% trackingStruct2.SD = trackingStruct2.SD.';
% trackingStruct2.d10 = trackingStruct2.d10.';
% trackingStruct2.d50 = trackingStruct2.d50.';
% trackingStruct2.d90 = trackingStruct2.d90.';
% trackingStruct2.ValidTracks = trackingStruct2.ValidTracks.';

%%
trackingStruct2.AvgDiffCoeff_ALL = trackingStruct2.AvgDiffCoeff_ALL.';
trackingStruct2.AvgDiffCoeff_FALSE = trackingStruct2.AvgDiffCoeff_FALSE.';
trackingStruct2.AvgDiffCoeff_TRUE = trackingStruct2.AvgDiffCoeff_TRUE.';
trackingStruct2.AvgLogIntensity_ALL = trackingStruct2.AvgLogIntensity_ALL.';
trackingStruct2.AvgLogIntensity_FALSE = trackingStruct2.AvgLogIntensity_FALSE.';
trackingStruct2.AvgLogIntensity_TRUE = trackingStruct2.AvgLogIntensity_TRUE.';
trackingStruct2.AvgTrackLength_ALL = trackingStruct2.AvgTrackLength_ALL.';
trackingStruct2.AvgTrackLength_FALSE = trackingStruct2.AvgTrackLength_FALSE.';
trackingStruct2.AvgTrackLength_TRUE = trackingStruct2.AvgTrackLength_TRUE.';
trackingStruct2.AvgSize_ALL = trackingStruct2.AvgSize_ALL.';
trackingStruct2.AvgSize_FALSE = trackingStruct2.AvgSize_FALSE.';
trackingStruct2.AvgSize_TRUE = trackingStruct2.AvgSize_TRUE.';
%% Assign Bead/Sample
for i = 1:length(trackingStruct2.Sample)
    if contains(trackingStruct2.Sample(i), "Beads") || contains(trackingStruct2.Sample(i), "QC") || contains(trackingStruct2.Sample(i), "100nm")
        trackingStruct2.SampleType(i,1) = "Beads";
    else
        trackingStruct2.SampleType(i,1) = "Sample";
    end
end

%% filter out fNTA/Scatter Runs
for i = 1:length(trackingStruct2.Sample)
    if contains(trackingStruct2.Sample(i), "fNTA") || contains(trackingStruct2.Sample(i), "Scatter") || contains(trackingStruct2.Sample(i), "nanoFCM") || contains(trackingStruct2.Sample(i), "fTNA")
        trackingStruct2.Sample(i,1) = missing; 
    end
end
%%
trackingTable = struct2table(trackingStruct2);
%remove rows containing missing data --> just rmmissing()
trackingTable_clean = rmmissing(trackingTable);
writetable(trackingTable_clean,'backgroundTracking-' + date + '.xlsx')