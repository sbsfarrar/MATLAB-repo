function [toBeProcessed, folder, answer, fileList] = fileSearch()
%fileSearch This function searches NAS for microscope images to be
%processed and checks them against processed files in the PD directory
%INPUT is the pwd that will be returned to when file search is complete
%OUTPUT is a nested structure of folder name and its accompanying file
%data
    toBeProcessed = struct; 
    homeEnvironment = pwd; %should be the Process Development Directory
    adsEnvironment = "";
    nasEnvironment = "Z:\PO SUSPENSION CULTURES";
    answer = questdlg('From NAS--> Select Lot Folder',...
        'Select Folder',...
        'SELECT', 'Access SharePoint folders','ADS Images','SELECT');
            switch answer
                case 'SELECT'
	                % Get the name of the folder that the user wants to use.
	                defaultFileName = fullfile(nasEnvironment, '*.*');
	                [folder] = uigetdir();
                    disp([folder ' was selected.'])
                case 'Access SharePoint folders'
	                cd(homeEnvironment)
                    defaultFileName = fullfile(pwd, '*.*');
 		            [folder] = uigetdir(defaultFileName, 'Select Image File Folder');
                    fileParse = sprintf('*%s*.tif',"4X");
                    fileList = dir(fullfile(folder,'**',fileParse));
                    return;
                case 'ADS Images'
                    %cd ..\..\ %'Process Development - P1-22006'\
                    %this directory shows all PD folders for you to choose,
                    %if 22006 is not present, then the folder is not synced
                    defaultFileName = fullfile(pwd, '*.*');
 		            [folder] = uigetdir(defaultFileName, 'Select Image File Folder');
                    fileParse = sprintf('.jpg');
                    fileList = dir(fullfile(folder,'**',fileParse));
                    return;
            end
    %%% Folder is selected, send to command window %%%
    %%% Search subdirectories until files are found %%%
    % filter out directory names with period indexes
    %A = dir(folder);
        %filter out period indexing
    %A = A(~startsWith({A.name},'.'));

    fileParse = sprintf('*%s*.tif',"4X");
    fileList = dir(fullfile(folder,'**',fileParse));

    %change back to Process Development Environment 
    %cd(homeEnvironment);

    %sort through fileList and create local directories of the sorted files
    for i = 1:numel(fileList)
        %check if processed file exists in the PD 
        %create folder name 
        vesselFileName = fileList(i).name;
        vesselStringSplit = strsplit(vesselFileName, ' ');
        vesselString = vesselStringSplit{1};
        vesselLetter = regexp(vesselString,'([A-Z]+)','match');
        folderDateString = []; 
        folderDateString = strsplit(fileList(i).folder,'\');
        folderDateString = folderDateString{end};
        %%%check if year is present, because people forget
        if contains(folderDateString,'2022') || contains(folderDateString,'2023')
            %miraculously all good! Do nothing
        else
            if contains(fileList(i).folder,'2022') %check upper directories for year 2022
                year = '2022';
            elseif contains(fileList(i).folder,'2023') %check upper directories for year 2023
                year = '2023';
            end
            folderDateString = [folderDateString year];
        end

        if isempty(vesselLetter) == 1 && numel(vesselString) == 4 %where zero is written instead of 'O'
            vesselString(1) = 'O';
            vesselLetter = regexp(vesselString,'([A-Z]+)','match'); %second try
            %fileList(i).name(1) = 'O';

            vesselLetter = vesselLetter{1};
        
            folderString = strsplit(folder,'\');
            folderString = folderString{end};
            folderString = strrep(folderString,'-','_');
            folderName = [folderString ' ' vesselLetter];
            folderFieldName = [folderString '_' vesselLetter];
    
            %folderDateString = strsplit(fileList(i).folder,'\');
            %folderDateString = folderDateString{end};
            %%% Search for typos in date, most commonly O and 0
            dateString = folderDateString(3:5);
            dateString(dateString == '0') = 'O';
            folderDateString(3:5) = dateString; 

            filePathString = [fileList(i).folder '\' vesselFileName];
            %continue; 
        elseif isempty(vesselLetter) == 1 && numel(vesselString) == 5 %where zero is written instead of 'O'
            vesselString(1:2) = 'OO';
            vesselLetter = regexp(vesselString,'([A-Z]+)','match'); %second try
            %fileList(i).name(1:2) = 'OO';

            vesselLetter = vesselLetter{1};
        
            folderString = strsplit(folder,'\');
            folderString = folderString{end};
            folderString = strrep(folderString,'-','_');
            folderName = [folderString ' ' vesselLetter];
            folderFieldName = [folderString '_' vesselLetter];
    
            %folderDateString = strsplit(fileList(i).folder,'\');
            %folderDateString = folderDateString{end};
            %%% Search for typos in date, most commonly O and 0
            dateString = folderDateString(3:5);
            dateString(dateString == '0') = 'O';
            folderDateString(3:5) = dateString;  
    
            filePathString = [fileList(i).folder '\' vesselFileName];
            %continue; 
        elseif string(vesselLetter) == "MM"
            if numel(vesselString) == 4 && ~contains(vesselString,'.')
                                                                            %capture specific case of MM19
                vesselString = [vesselString(1:2) '1' vesselString(3:4)];
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
    
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif numel(vesselString) == 6 && ~contains(vesselString,'.') 
                                                                                %capture specific case of MM1221
                vesselString(end) = []; 
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
    
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif numel(vesselString) == 6 && contains(vesselString,'.')
                %capture specific case of MM114., MM14.5, and MM12.1
                additionalSplit = strsplit(vesselString,'.');
                additionalSplit = additionalSplit{1};
                if numel(additionalSplit) == 4 && string(folderDateString) == "07DEC2022"
                    additionalSplit(4:5) = '24';
                    vesselString = additionalSplit; 
                    %vesselStringSplit{1} = additionalSplit; 
                    %newVesselString = strjoin(vesselStringSplit,' ');
                elseif numel(additionalSplit) == 4 && string(folderDateString) == "08DEC2022"
                    additionalSplit(4:5) = '25';
                    vesselString = additionalSplit;
                    %vesselStringSplit{1} = additionalSplit; 
                    %newVesselString = strjoin(vesselStringSplit,' ');
                else
                    vesselString(end) = []; 
                    vesselString(4) = '2';
                    %vesselStringSplit{1} = vesselString; 
                    %vesselStringSplit(2) = []; %remove extraneous integer in file name
                    %newVesselString = strjoin(vesselStringSplit,' ');
                    %fileList(i).name = newVesselString;
                end

                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif numel(vesselString) == 7 && contains(vesselString,'114')
                                                                                %capture specific case of MM114.#
                vesselString(6:end) = []; 
                vesselString(4) = '2';
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
    
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif numel(vesselString) == 7 ...
                    && numel(vesselLetter) > 1              %capture specific case of MM1.14X and MM1.34X
                                                            %supposed to be MM123.1 & MM123.3, change to MM123                                                      
                vesselString(4:5) = '23';
                vesselString(6:end) = [];
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
    
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif numel(vesselString) == 4 && contains(vesselString,'.')             
                                                            %capture specific case of MM.1 and MM.3
                                                            %supposed to be MM123.1 & MM123.3, change to MM123                                                      
                vesselString(3:5) = '123';
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;
            elseif string(folderDateString) == "14NOV2022"
                                                                                %capture specific case of wrong vessel ID
                vesselString = 'MM122';
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
    
                vesselLetter = vesselLetter{1};
                
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;  
            else 
                vesselLetter = vesselLetter{1};
        
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' fileList(i).name];
            end
        elseif string(vesselLetter) == "VB" %case where VB is supposed to be B
                vesselLetter = 'B';
                vesselString(1) = [];
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;  
        elseif string(vesselLetter) == "BB" %case where BB is supposed to be B
                vesselLetter = 'B';
                vesselString(1) = [];
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;  
         elseif string(vesselLetter) == "XC" %case where XC is supposed to be X
                vesselLetter = 'X';
                vesselString(2) = [];
                %vesselStringSplit{1} = vesselString; 
                %newVesselString = strjoin(vesselStringSplit,' ');
                %fileList(i).name = newVesselString;
                folderString = strsplit(folder,'\');
                folderString = folderString{end};
                folderString = strrep(folderString,'-','_');
                folderName = [folderString ' ' vesselLetter];
                folderFieldName = [folderString '_' vesselLetter];
        
                %folderDateString = strsplit(fileList(i).folder,'\');
                %folderDateString = folderDateString{end};
                %%% Search for typos in date, most commonly O and 0
                dateString = folderDateString(3:5);
                dateString(dateString == '0') = 'O';
                folderDateString(3:5) = dateString;  
        
                filePathString = [fileList(i).folder '\' vesselFileName];
                %continue;  
        else 
            try %attempt to catch lower-case letters
                vesselLetter = vesselLetter{1};
            catch
                vesselString = upper(vesselString);
                vesselLetter = regexp(vesselString,'([A-Z]+)','match');
                vesselLetter = vesselLetter{1};
            end

            folderString = strsplit(folder,'\');
            folderString = folderString{end};
            folderString = strrep(folderString,'-','_');
            folderName = [folderString ' ' vesselLetter];
            folderFieldName = [folderString '_' vesselLetter];
    
            %folderDateString = strsplit(fileList(i).folder,'\');
            %folderDateString = folderDateString{end};
            %%% Search for typos in date, most commonly O and 0
            dateString = folderDateString(3:5);
            dateString(dateString == '0') = 'O';
            folderDateString(3:5) = dateString;  
    
            filePathString = [fileList(i).folder '\' fileList(i).name];
        end
        disp(folderFieldName)
        if isfield(toBeProcessed,folderFieldName) == 0
          toBeProcessed.(folderFieldName) = [];
          toBeProcessed.(folderFieldName).FileName = [];
          toBeProcessed.(folderFieldName).AdjustedFileName = [];
          toBeProcessed.(folderFieldName).ProcessedFileName = [];
          toBeProcessed.(folderFieldName).FileDate = [];
          toBeProcessed.(folderFieldName).FilePath = [];
        end

        toBeProcessed.(folderFieldName).FileName = [toBeProcessed.(folderFieldName).FileName; convertCharsToStrings(fileList(i).name)];
        toBeProcessed.(folderFieldName).FileDate = [toBeProcessed.(folderFieldName).FileDate; convertCharsToStrings(folderDateString)];
        toBeProcessed.(folderFieldName).FilePath = [toBeProcessed.(folderFieldName).FilePath; convertCharsToStrings(filePathString)];

        %create processed and adjusted file name --> used for saving processing files and
        %checking if existing files are already processed.

        % IF Vessels contain '.' annotation, remove '.' and treat as
        % regular ID
        if contains(vesselString,'.') == 1
            vesselSplit = strsplit(vesselString,'.');
            vesselString = vesselSplit{1};
        end

        adjustedFileName = [vesselString ' ' folderDateString ' 4X']; %no extension
        processedFileName = [vesselString '-' folderDateString]; %no extension
        if i > 1 && isempty(toBeProcessed.(folderFieldName).AdjustedFileName) == 0
            strAdj = convertStringsToChars(toBeProcessed.(folderFieldName).AdjustedFileName);
            strPro = convertStringsToChars(toBeProcessed.(folderFieldName).ProcessedFileName);
            adjSearch = strfind(strAdj, adjustedFileName,"ForceCellOutput",true);
            proSearch = strfind(strPro, processedFileName,"ForceCellOutput",true);
            adjSearchEmpty = cellfun(@isempty,adjSearch);
            proSearchEmpty = cellfun(@isempty,proSearch);
            adjSearchCount = find(adjSearchEmpty==0);
            proSearchCount = find(proSearchEmpty==0);
            if numel(adjSearchCount) >= 1
                x = num2str(numel(adjSearchCount));
                adjustedFileName = [vesselString ' ' folderDateString ' 4X (' x ').tif'];
            else
                adjustedFileName = [vesselString ' ' folderDateString ' 4X.tif'];
            end
            if numel(proSearchCount) >= 1
                y = num2str(numel(proSearchCount));
                processedFileName = [vesselString '-' folderDateString ' (' y ').mat'];
            else
                processedFileName = [vesselString '-' folderDateString '.mat'];
            end
        else 
            adjustedFileName = [vesselString ' ' folderDateString ' 4X.tif'];
            processedFileName = [vesselString '-' folderDateString '.mat'];
        end
        
        toBeProcessed.(folderFieldName).AdjustedFileName = [toBeProcessed.(folderFieldName).AdjustedFileName; convertCharsToStrings(adjustedFileName)];
        toBeProcessed.(folderFieldName).ProcessedFileName = [toBeProcessed.(folderFieldName).ProcessedFileName; convertCharsToStrings(processedFileName)];

    end

    %%%% LOOK FOR TYPOS BEFORE PROCESSING %%%%
    %%% separate script will be generated to better organize errors %%%





end