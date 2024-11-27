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
        date_str = string(regexp(runString,'[0-9]+-[0-9]+-[0-9]+','match'));
        sample_summ = strrep(images(i).name,"-background.bmp","");
        sample_summ = strrep(sample_summ,regexp(sample_summ,"[0-9]+-[0-9]+-[0-9]+ [0-9]+-[0-9]+-[0-9]+","match"),"");
        sample_summ = strtrim(sample_summ);
        imagePath = images(i).folder + "\" + images(i).name;
        trackingStruct.Sample(i,1) = sample_summ; 
        trackingStruct.Date(i,1) = date_str; 
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
nta_path = "backgroundTracking_full.xlsx";
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
%% if re-initializing
trackingStruct_filt = trackingStruct; 
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
            opts = detectImportOptions(summaryPath);
            opts = setvartype(opts, 'string');
            opts.DataLines = [1, 55; 58, 59; 68, 87]; %grabs everything but bin data, removing ppf warning
            % summary_table = readtable(summaryPath,opts);
            % Vars = table2array(summary_table(:,1)).';
            % Vals = table2array(summary_table(:,2)).';
            % adj_summary_table = array2table(Vals, 'VariableNames', Vars);
            try
                summary_table = readtable(summaryPath,opts);
                Vars = table2array(summary_table(:,1)).';
                var_names = summary_table.Properties.VariableNames;
                if length(var_names) > 2
                    Vals = table2array(summary_table(:,2:length(var_names)-1)).';
                else
                    Vals = table2array(summary_table(:,2)).';
                end
                adj_summary_table = array2table(Vals, 'VariableNames', Vars);
                trackingStruct_filt.Concentration(i,1) = double(adj_summary_table.("Concentration (Particles / ml)")(1));
                trackingStruct_filt.PPF(i,1) = double(adj_summary_table.("Particles per frame")(1));
                trackingStruct_filt.CPF(i,1) = double(adj_summary_table.("Centres per frame")(1));
                trackingStruct_filt.CompletedTracks(i,1) = double(adj_summary_table.("Completed tracks")(1));
                trackingStruct_filt.xDrift(i,1) = double(adj_summary_table.("X-Drift (pix/frame)")(1));
                trackingStruct_filt.yDrift(i,1) = double(adj_summary_table.("Y-Drift (pix/frame)")(1));
                trackingStruct_filt.Mean(i,1) = double(adj_summary_table.Mean(1));
                trackingStruct_filt.Mode(i,1) = double(adj_summary_table.Mode(1));
                trackingStruct_filt.SD(i,1) = double(adj_summary_table.SD(1));
                trackingStruct_filt.d10(i,1) = double(adj_summary_table.D10(1));
                trackingStruct_filt.d50(i,1) = double(adj_summary_table.D50(1));
                trackingStruct_filt.d90(i,1) = double(adj_summary_table.D90(1));
                trackingStruct_filt.ValidTracks(i,1) = double(adj_summary_table.("Valid Tracks")(1));

                trackingStruct_filt.SoftwareVersion(i,1) = adj_summary_table.("Software Version")(1);
                trackingStruct_filt.Diluent(i,1) = adj_summary_table.Diluent(1);
                if length(var_names) > 2
                    agg_string = "";
                    for k = 1:length(adj_summary_table.Remarks)
                        agg_string = agg_string + " " + adj_summary_table.Remarks(k);
                    end
                        remarks_split = strsplit(agg_string);
                        idx_focus = find(remarks_split == "Focus");
                else
                    remarks_split = strsplit(adj_summary_table.Remarks);
                    idx_focus = find(remarks_split == "Focus");
                end
                if isempty(idx_focus)
                    trackingStruct_filt.Focus(i,1) = 999999999;
                else
                    trackingStruct_filt.Focus(i,1) = double(regexp(remarks_split(idx_focus+1),'-?[\d]+','match'));
                end
                trackingStruct_filt.Temperature(i,1) = adj_summary_table.("Temperature/C")(1);
                trackingStruct_filt.Viscosity(i,1) = adj_summary_table.("Viscosity/cP")(1);
                trackingStruct_filt.CameraType(i,1) = adj_summary_table.("Camera Type")(1);
                trackingStruct_filt.LaserType(i,1) = adj_summary_table.("Laser Type")(1);
                trackingStruct_filt.CameraLevel(i,1) = adj_summary_table.("Camera Level")(1);
                trackingStruct_filt.SliderShutter(i,1) = adj_summary_table.("Slider Shutter")(1);
                trackingStruct_filt.SliderGain(i,1) = adj_summary_table.("Slider Gain")(1);
                trackingStruct_filt.Shutter(i,1) = adj_summary_table.("Shutter/ms")(1);
                trackingStruct_filt.CameraHistUL(i,1) = adj_summary_table.("Camera Histogram Upper Limit")(1);
                trackingStruct_filt.CameraHistLL(i,1) = adj_summary_table.("Camera Histogram Lower Limit")(1);
                trackingStruct_filt.FrameRate(i,1) = adj_summary_table.("Frame rate/fps")(1);
                trackingStruct_filt.PumpSpeed(i,1) = adj_summary_table.("Syringe Pump Speed/AU")(1);
                trackingStruct_filt.DetectionThreshold(i,1) = adj_summary_table.("Detection Threshold")(1);
                trackingStruct_filt.MaxJumpDist(i,1) = adj_summary_table.("Max Jump Distance")(1);
                trackingStruct_filt.TotalFrames(i,1) = adj_summary_table.("Total frames analysed")(1);
                trackingStruct_filt.DilutionFactor(i,1) = double(adj_summary_table.("Dilution factor (concentrations adjusted for this factor)")(1));
                
                trackingStruct_filt.Concentration_Warning(i,1) = adj_summary_table.("Validity of concentration measurement")(1);
                trackingStruct_filt.NoiseLevel_Warning(i,1) = adj_summary_table.("Noise level")(1);
                trackingStruct_filt.Brightness_Warning(i,1) = adj_summary_table.("Brightness")(1);
                trackingStruct_filt.Vibration_Warning(i,1) = adj_summary_table.("Vibration detected")(1);
                trackingStruct_filt.AnalysisMethod(i,1) = adj_summary_table.("Analysis Method")(1);
                trackingStruct_filt.Weighting(i,1) = adj_summary_table.("Weighting")(1);

                disp("Summary Loaded")
            catch
                trackingStruct_filt.Concentration(i,1) = missing;
                trackingStruct_filt.PPF(i,1) = missing;
                trackingStruct_filt.CPF(i,1) = missing;
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

                trackingStruct_filt.SoftwareVersion(i,1) = missing;
                trackingStruct_filt.Diluent(i,1) = missing;
                trackingStruct_filt.Focus(i,1) = missing;
                trackingStruct_filt.Temperature(i,1) = missing;
                trackingStruct_filt.Viscosity(i,1) = missing;
                trackingStruct_filt.CameraType(i,1) = missing;
                trackingStruct_filt.LaserType(i,1) = missing;
                trackingStruct_filt.CameraLevel(i,1) = missing;
                trackingStruct_filt.SliderShutter(i,1) = missing;
                trackingStruct_filt.SliderGain(i,1) = missing;
                trackingStruct_filt.Shutter(i,1) = missing;
                trackingStruct_filt.CameraHistUL(i,1) = missing;
                trackingStruct_filt.CameraHistLL(i,1) = missing;
                trackingStruct_filt.FrameRate(i,1) = missing;
                trackingStruct_filt.PumpSpeed(i,1) = missing;
                trackingStruct_filt.DetectionThreshold(i,1) = missing;
                trackingStruct_filt.MaxJumpDist(i,1) = missing;
                trackingStruct_filt.TotalFrames(i,1) = missing;
                trackingStruct_filt.DilutionFactor(i,1) = missing;
                
                trackingStruct_filt.Concentration_Warning(i,1) = missing;
                trackingStruct_filt.NoiseLevel_Warning(i,1) = missing;
                trackingStruct_filt.Brightness_Warning(i,1) = missing;
                trackingStruct_filt.Vibration_Warning(i,1) = missing;
                trackingStruct_filt.AnalysisMethod(i,1) = missing;
                trackingStruct_filt.Weighting(i,1) = missing;
                disp("skipping summary")
            end
            break;
        else
            trackingStruct_filt.Concentration(i,1) = 999999999;
            trackingStruct_filt.PPF(i,1) = 999999999;
            trackingStruct_filt.CPF(i,1) = 999999999;
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

            trackingStruct_filt.SoftwareVersion(i,1) = "";
            trackingStruct_filt.Diluent(i,1) = "";
            trackingStruct_filt.Focus(i,1) = 999999999;
            trackingStruct_filt.Temperature(i,1) = "";
            trackingStruct_filt.Viscosity(i,1) = "";
            trackingStruct_filt.CameraType(i,1) = "";
            trackingStruct_filt.LaserType(i,1) = "";
            trackingStruct_filt.CameraLevel(i,1) = "";
            trackingStruct_filt.SliderShutter(i,1) = "";
            trackingStruct_filt.SliderGain(i,1) = "";
            trackingStruct_filt.Shutter(i,1) = "";
            trackingStruct_filt.CameraHistUL(i,1) = "";
            trackingStruct_filt.CameraHistLL(i,1) = "";
            trackingStruct_filt.FrameRate(i,1) = "";
            trackingStruct_filt.PumpSpeed(i,1) = "";
            trackingStruct_filt.DetectionThreshold(i,1) = "";
            trackingStruct_filt.MaxJumpDist(i,1) = "";
            trackingStruct_filt.TotalFrames(i,1) = "";
            trackingStruct_filt.DilutionFactor(i,1) = 999999999;
            
            trackingStruct_filt.Concentration_Warning(i,1) = "";
            trackingStruct_filt.NoiseLevel_Warning(i,1) = "";
            trackingStruct_filt.Brightness_Warning(i,1) = "";
            trackingStruct_filt.Vibration_Warning(i,1) = "";
            trackingStruct_filt.AnalysisMethod(i,1) = "";
            trackingStruct_filt.Weighting(i,1) = "";
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %UNDER DEVELOPMENT -- EXPERIMENT REPORTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runQC = dir(dirPathQC + "\*\*-ExperimentSummary.csv*");
runPD = dir(dirPathPD + "\*\*-ExperimentSummary.csv*");
runArch = dir(dirPathArch + "\*\*-ExperimentSummary.csv*");
run_full = [runQC; runPD; runArch];

%remove 'Combined experiment'
A = struct2cell(run_full).';
B = cellfun('isempty',strfind(A(:,1), "Combined experiment")); %logical index, yay!
run_full = run_full(B); 
C = struct2cell(run_full).';
D = cellfun('isempty',strfind(C(:,1), "Combined Experiment")); %logical index, yay!
run_full = run_full(D);

%initialize new stucture variables
trackingStruct_filt.RunName = strings(length(trackingStruct_filt.Sample), 1);
trackingStruct_filt.RunConcentration_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunConcentration_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunConcentration_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunConcentration_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunSize_Mean_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mean_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mean_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mean_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunSize_Mode_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mode_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mode_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_Mode_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunSize_d10_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d10_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d10_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d10_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunSize_d50_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d50_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d50_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d50_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunSize_d90_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d90_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d90_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunSize_d90_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

trackingStruct_filt.RunValidTracks_Mean = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunValidTracks_StDev = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunValidTracks_StErr = ones(length(trackingStruct_filt.Sample), 1)*999999999;
trackingStruct_filt.RunValidTracks_PercentCV = ones(length(trackingStruct_filt.Sample), 1)*999999999;

for j=1:length(run_full)
    disp(j)
    %disp(length(summary) - 1);
    runPath = run_full(j).folder + "\" + run_full(j).name;
    sample_run_split = strsplit(run_full(j).name);
    sample_run = strtrim(strjoin(sample_run_split(1:end-1), " ")); %remove individual appended items (keep date)
    %sample_run = strtrim(sample_run);
    %find index of sample_run in trackingStruct_filt
    idx = ~cellfun('isempty',strfind(trackingStruct_filt.Path, sample_run)); %logical index, yay!
    if any(idx) && any(trackingStruct_filt.RunName(idx) == "")
        disp('found!')
        opts = detectImportOptions(runPath);
        opts = setvartype(opts, 'string');
        opts.DataLines = [1, 55; 58, 59; 68, 87];
        run_table = readtable(runPath,opts);
    
        % reformat table
        Vars = table2array(run_table(:,1)).';
        var_names = run_table.Properties.VariableNames;
        if length(var_names) > 2
            Vals = table2array(run_table(:,2:length(var_names)-1)).';
        else
            Vals = table2array(run_table(:,2)).';
        end

        adj_run_table = array2table(Vals, 'VariableNames', Vars);
        trim_max = length(rmmissing(adj_run_table.("Camera Type")));
        trimmed_run_table = adj_run_table(1:trim_max,:); 
 
        try
            %populate structure, using logical index variable
            trackingStruct_filt.RunName(idx) = sample_run;
        
            trackingStruct_filt.RunConcentration_Mean(idx) = mean(double(trimmed_run_table.("Concentration (Particles / ml)")));
            trackingStruct_filt.RunConcentration_StDev(idx) = std(double(trimmed_run_table.("Concentration (Particles / ml)")));
            trackingStruct_filt.RunConcentration_StErr(idx) = std(double(trimmed_run_table.("Concentration (Particles / ml)")))/...
                sqrt(length(double(trimmed_run_table.("Concentration (Particles / ml)"))));
            trackingStruct_filt.RunConcentration_PercentCV(idx) = (std(double(trimmed_run_table.("Concentration (Particles / ml)")))/...
                mean(double(trimmed_run_table.("Concentration (Particles / ml)"))))*100;
            
            trackingStruct_filt.RunSize_Mean_Mean(idx) = mean(double(trimmed_run_table.Mean));
            trackingStruct_filt.RunSize_Mean_StDev(idx) = std(double(trimmed_run_table.Mean));
            trackingStruct_filt.RunSize_Mean_StErr(idx) = std(double(trimmed_run_table.Mean))/...
                sqrt(length(double(trimmed_run_table.Mean)));
            trackingStruct_filt.RunSize_Mean_PercentCV(idx) = (std(double(trimmed_run_table.Mean))/...
                mean(double(trimmed_run_table.Mean)))*100;
            
            trackingStruct_filt.RunSize_Mode_Mean(idx) = mean(double(trimmed_run_table.Mode));
            trackingStruct_filt.RunSize_Mode_StDev(idx) = std(double(trimmed_run_table.Mode));
            trackingStruct_filt.RunSize_Mode_StErr(idx) = std(double(trimmed_run_table.Mode))/...
                sqrt(length(double(trimmed_run_table.Mode)));
            trackingStruct_filt.RunSize_Mode_PercentCV(idx) = (std(double(trimmed_run_table.Mode))/...
                mean(double(trimmed_run_table.Mode)))*100;
            
            trackingStruct_filt.RunSize_d10_Mean(idx) = mean(double(trimmed_run_table.D10));
            trackingStruct_filt.RunSize_d10_StDev(idx) = std(double(trimmed_run_table.D10));
            trackingStruct_filt.RunSize_d10_StErr(idx) = std(double(trimmed_run_table.D10))/...
                sqrt(length(double(trimmed_run_table.D10)));
            trackingStruct_filt.RunSize_d10_PercentCV(idx) = (std(double(trimmed_run_table.D10))/...
                mean(double(trimmed_run_table.D10)))*100;
            
            trackingStruct_filt.RunSize_d50_Mean(idx) = mean(double(trimmed_run_table.D50));
            trackingStruct_filt.RunSize_d50_StDev(idx) = std(double(trimmed_run_table.D50));
            trackingStruct_filt.RunSize_d50_StErr(idx) = std(double(trimmed_run_table.D50))/...
                sqrt(length(double(trimmed_run_table.D50)));
            trackingStruct_filt.RunSize_d50_PercentCV(idx) = (std(double(trimmed_run_table.D50))/...
                mean(double(trimmed_run_table.D50)))*100;
            
            trackingStruct_filt.RunSize_d90_Mean(idx) = mean(double(trimmed_run_table.D90));
            trackingStruct_filt.RunSize_d90_StDev(idx) = std(double(trimmed_run_table.D90));
            trackingStruct_filt.RunSize_d90_StErr(idx) = std(double(trimmed_run_table.D90))/...
                sqrt(length(double(trimmed_run_table.D90)));
            trackingStruct_filt.RunSize_d90_PercentCV(idx) = (std(double(trimmed_run_table.D90))/...
                mean(double(trimmed_run_table.D90)))*100;
            
            trackingStruct_filt.RunValidTracks_Mean(idx) = mean(double(trimmed_run_table.("Valid Tracks")));
            trackingStruct_filt.RunValidTracks_StDev(idx) = std(double(trimmed_run_table.("Valid Tracks")));
            trackingStruct_filt.RunValidTracks_StErr(idx) = std(double(trimmed_run_table.("Valid Tracks")))/...
                sqrt(length(double(trimmed_run_table.("Valid Tracks"))));
            trackingStruct_filt.RunValidTracks_PercentCV(idx) = (std(double(trimmed_run_table.("Valid Tracks")))/...
                mean(double(trimmed_run_table.("Valid Tracks"))))*100;
        catch
            names = fieldnames(trackingStruct_filt);
            for k = 1:length(names)
                if contains(names{k}, "Run")
                    if class(trackingStruct_filt.(names{k})) == "string"
                        trackingStruct_filt.(names{k}) = "";
                    else
                        trackingStruct_filt.(names{k}) = 999999999; 
                    end
                end
            end
        end
    end
end

%% for adjusting original nta particle data

trackingStruct2 = table2struct(nta_particle_data, "ToScalar", true);
trackingStruct2 = rmfield(trackingStruct2, "SampleType");
names = fieldnames(trackingStruct2);
for i = 1:length(names)
    if class(trackingStruct2.(names{i})) == "cell"
        trackingStruct2.(names{i}) = string(trackingStruct2.(names{i}));
    end
end
%% FINAL CLEAN AND EXPORT
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

% trackingStruct2.AvgDiffCoeff_ALL = trackingStruct2.AvgDiffCoeff_ALL.';
% trackingStruct2.AvgDiffCoeff_FALSE = trackingStruct2.AvgDiffCoeff_FALSE.';
% trackingStruct2.AvgDiffCoeff_TRUE = trackingStruct2.AvgDiffCoeff_TRUE.';
% trackingStruct2.AvgLogIntensity_ALL = trackingStruct2.AvgLogIntensity_ALL.';
% trackingStruct2.AvgLogIntensity_FALSE = trackingStruct2.AvgLogIntensity_FALSE.';
% trackingStruct2.AvgLogIntensity_TRUE = trackingStruct2.AvgLogIntensity_TRUE.';
% trackingStruct2.AvgTrackLength_ALL = trackingStruct2.AvgTrackLength_ALL.';
% trackingStruct2.AvgTrackLength_FALSE = trackingStruct2.AvgTrackLength_FALSE.';
% trackingStruct2.AvgTrackLength_TRUE = trackingStruct2.AvgTrackLength_TRUE.';
% trackingStruct2.AvgSize_ALL = trackingStruct2.AvgSize_ALL.';
% trackingStruct2.AvgSize_FALSE = trackingStruct2.AvgSize_FALSE.';
% trackingStruct2.AvgSize_TRUE = trackingStruct2.AvgSize_TRUE.';
%% Assign Bead/Sample
for i = 1:length(trackingStruct2.Sample)
    if contains(trackingStruct2.Sample(i), "Beads") || contains(trackingStruct2.Sample(i), "QC") || ...
            contains(trackingStruct2.Sample(i), "100nm")
        trackingStruct2.SampleType(i,1) = "Beads";
    elseif contains(trackingStruct2.Sample(i), "Vehicle") || contains(trackingStruct2.Sample(i), "Placebo")
        trackingStruct2.SampleType(i,1) = "Vehicle";
    elseif contains(trackingStruct2.Sample(i), "fNTA") || contains(trackingStruct2.Sample(i), "Scatter") || ...
                contains(trackingStruct2.Sample(i), "nanoFCM") || contains(trackingStruct2.Sample(i), "fTNA") || ...
                contains(trackingStruct2.Sample(i), "Alexa488") || contains(trackingStruct2.Sample(i), "CD73") || ...
                contains(trackingStruct2.Sample(i), "ExoGlow") || contains(trackingStruct2.Sample(i), "CFSE") || ...
                contains(trackingStruct2.Sample(i), "MemGlow488")
            
            trackingStruct2.SampleType(i,1) = "fNTA"; 
    else
        trackingStruct2.SampleType(i,1) = "Sample";
    end
end
%% Assign Sample Class
for i = 1:length(trackingStruct2.Sample)
    if trackingStruct2.SampleType(i) == "Beads"
        trackingStruct2.SampleClass(i,1) = "Beads";
    elseif trackingStruct2.SampleType(i) == "Vehicle"
        trackingStruct2.SampleClass(i,1) = "Vehicle";
    elseif trackingStruct2.SampleType(i) == "fNTA"
        trackingStruct2.SampleClass(i,1) = "fNTA";
    elseif trackingStruct2.SampleType(i) == "Sample"
        if contains(trackingStruct2.Sample(i), "IN1")
            trackingStruct2.SampleClass(i,1) = "IN1";
        elseif (contains(trackingStruct2.Sample(i), "IN2") && contains(trackingStruct2.Sample(i), "pre")) || ...
                contains(trackingStruct2.Sample(i), "DI1")
            trackingStruct2.SampleClass(i,1) = "IN2pre";
        elseif contains(trackingStruct2.Sample(i), "IN2") && contains(trackingStruct2.Sample(i), "post")
            trackingStruct2.SampleClass(i,1) = "IN2post";
        elseif contains(trackingStruct2.Sample(i), "DI2")
            trackingStruct2.SampleClass(i,1) = "DI2";
        elseif contains(trackingStruct2.Sample(i), "IN3")
            trackingStruct2.SampleClass(i,1) = "IN3";
        elseif contains(trackingStruct2.Sample(i), "DS3")
            trackingStruct2.SampleClass(i,1) = "DS3";
        elseif contains(trackingStruct2.Sample(i), "DP")
            trackingStruct2.SampleClass(i,1) = "DP";
        else
            trackingStruct2.SampleClass(i,1) = "Other";
        end
    end
end
%% filter out mean value errors (zeros, 999999999, and NaN)
names = fieldnames(trackingStruct2);
for i=1:length(names)
    % filter out fNTA/Scatter Runs
    if class(trackingStruct2.(names{i})) == "string"
        %names{i} == "Date" || names{i} == "Operator" || names{i} == "Path" || names{i} == "SampleType"
        %do nothing
        for j = 1:length(trackingStruct2.(names{i}))
            if trackingStruct2.(names{i})(j) == ""
                trackingStruct2.(names{i})(j) = missing; 
            end
        end
    elseif class(trackingStruct2.(names{i})) == "double" %some focus settings can be zero 
        if names{i} == "Focus" || names{i} == "PPF" || names{i} == "CPF" || ...
                names{i} == "xDrift" || names{i} == "yDrift" || contains(names{i}, "Run")
            for j = 1:length(trackingStruct2.(names{i}))
                if trackingStruct2.(names{i})(j) == 999999999
                    trackingStruct2.(names{i})(j) = missing; %will make them NaN for doubles
                end
            end
        else
            for j = 1:length(trackingStruct2.(names{i}))
                if trackingStruct2.(names{i})(j) == 0 || trackingStruct2.(names{i})(j) == 999999999
                    trackingStruct2.(names{i})(j) = missing; %will make them NaN for doubles
                end
            end
        end
    end
end

%% Export tables for JMP

trackingTable = struct2table(trackingStruct2);
%remove rows containing missing data --> just rmmissing()
%trackingTable_clean = rmmissing(trackingTable);
%trackingTable_full = [nta_particle_data; trackingTable_clean];
trackingTable_full = [nta_particle_data; trackingTable];
writetable(trackingTable,'backgroundTracking-' + string(datetime(("today"))) + '.xlsx')
%writetable(trackingTable_clean,'backgroundTracking.xlsx')
writetable(trackingTable_full,'backgroundTracking_full.xlsx') 


%%% Create backgroundTracking Suffix to Append
% suffixPath = "C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\MATLAB-repo\Nanoparticle Tracking Analysis";
% filesNTA = dir(suffixPath + "\*backgroundTracking*.xlsx");

