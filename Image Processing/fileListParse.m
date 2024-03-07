function [toBeProcessed, fileList, overwriteIndex] = fileListParse(fileList, toBeProcessed)
%fileListParse take the generated fileList from NAS along with the output
%structure for detection of files that have already been processed. 
%user will be prompted to either eliminate repeated files from the
%structure or overwrite them
%OUTPUT is the toBeProcessed structure and fileList array either undisturbed 
% (lame) or trimmed to include only novel files. 
% overwriteIndex will be a structure passed out to reference when specific
% files and their respective locations need to be replaced.

%first, perform a general directory and subdirectory search of all files in
%the root directory 
rootdir = pwd; %home environment (when published, this should be the process development folder)
               %or, for ADS images, the folder of the operator one
               %directory forward in the ADS process development folder 
processedFileList = dir(fullfile(rootdir,'**\*.mat')); %the processed sphere data is saved with 
               % an accompanying .mat file so we don't have to read and
               % reference an excel file each time. 
if numel(processedFileList) == 0
    return; %leaves function since there is no need to take action on fileList
end
%processedFileList = processedFileList(~[fileList.isdir]); %remove folders from the list
 %commented out the above statement because .mat* will automatically find
 %the files in question. 
numStr = num2str(numel(processedFileList)); 
msg = ['It appears there are ' numStr ' processed files present. Please select an option:'];
answer = questdlg(msg,...
    'Processed Files Present',...
    'Skip All Duplicate Files', 'Overwrite Files in Designated Folder(s)','Skip All Duplicate Files'); %Overwrite All no longer used
answerCase = 1;
overwriteIndex = 0; 
while answerCase == 1
    switch answer
        case 'Skip All Duplicate Files'
            %Now, we compare the processedFileList (the names fetched are adjusted
            %file names) to the toBeProcessed input structure
            fieldName = fieldnames(toBeProcessed);
            %repeatCount = 0;                
            for i = 1:numel(fieldName)
                subStructure = toBeProcessed.(fieldName{i});
                tableSubStructure = struct2table(subStructure);
                nonScalarSubStructure = table2struct(tableSubStructure);
                for j = 1 : numel(processedFileList)
                    processedFileListName = string(processedFileList(j).name);
                    %delete in fileList
                    repeatIndex = find(subStructure.ProcessedFileName == processedFileListName,1);
                    if isempty(repeatIndex)
                        %no file found, do nothing
                    else
                        repeatFileName = subStructure.FileName(repeatIndex); %string
                        repeatFolder = subStructure.FilePath(repeatIndex);
                        repeatFolder = strsplit(repeatFolder, '\');
                        repeatFolder = repeatFolder(1:end-1); %trim file name
                        repeatFolder = string(strjoin(repeatFolder,'\')); %piece back together for folder path
                        tempName = num2cell(string({fileList.name}));
                        tempFolder = num2cell(string({fileList.folder}));
                        [fileList.folder] = tempFolder{:};
                        [fileList.name] = tempName{:}; 
                        fileList([fileList.name] == repeatFileName & [fileList.folder] == repeatFolder) = [];
                    end
                    %delete in toBeProcessed
                    if isempty(nonScalarSubStructure)
                        break;
                    else
                        nonScalarSubStructure([nonScalarSubStructure.ProcessedFileName] == processedFileListName) = []; 
                    end
                end
                tableSubStructure = struct2table(nonScalarSubStructure);
                scalarSubStructure = table2struct(tableSubStructure,"ToScalar",true);
                if isempty(scalarSubStructure.FileName)
                    toBeProcessed = rmfield(toBeProcessed,fieldName{i});
                else
                    toBeProcessed.(fieldName{i}) = scalarSubStructure;
                end
                %empty nonScalarSubStructure
                nonScalarSubStructure = []; 
            end
            answerCase = 0; 
        case 'Overwrite Files in Designated Folder(s)'
            %prompt list dialogue to skip files not included in
            %selected list, return index of rowns to be replaced in
            %individual excel file (address aggregated data as well)
            folderPath = [rootdir '\Microscope Images\'];
            folderDirectory = dir(folderPath);
            folderDirectory = folderDirectory(~startsWith({folderDirectory.name},'.'));
            folderNames = {folderDirectory.name}';
            [indx,tf] = listdlg('PromptString',{'Select the Folder whose contents you ',...
                'wish to overwrite.',...
                'Only one folder can be selected at a time.',''},...
                'SelectionMode','single','ListString',folderNames); %If the user clicks Cancel, presses Esc, or clicks the 
                                           % close button (X) in the dialog box title bar, then the tf return value is 0. 
            if tf == 0
                uiwait(msgbox({'No selection was recorded.',...
                    'Returning to default: All files skipped.'}));
                answerCase = 1;
                answer = 'Skip All Duplicate Files';
                continue; 
            end
            newFolderPath = [folderPath folderNames{indx}];
            newFolderPathFiles = [newFolderPath '\Individual Image Data\'];
            fileParse = sprintf('*.mat');
            fileDirectory = dir(fullfile(newFolderPathFiles,'**',fileParse));
            fileNames = {fileDirectory.name}';
            %fileFolders = {fileDirectory.folder}';
            [indx2,tf2] = listdlg('PromptString',{'Select the files you wish ',...
                'to reprocess.'},...
                'ListString',fileNames); %If the user clicks Cancel, presses Esc, or clicks the 
            if tf2 == 0
                uiwait(msgbox({'No selection was recorded.',...
                    'Returning to original prompt.'}));
                answerCase = 1;
                continue; 
            else
                overwriteIndex = fileDirectory(indx2); %files to remove from processedFileList
                for i = 1: numel(overwriteIndex)
                    tempNameOverWrite = num2cell(string({overwriteIndex.name}));
                    tempNameProcessed = num2cell(string({processedFileList.name}));
                    [overwriteIndex.name] = tempNameOverWrite{:};
                    [processedFileList.name] = tempNameProcessed{:}; 
                    processedFileList([processedFileList.name] == overwriteIndex(i).name) = [];
                end
                answerCase = 1; 
                answer = 'Skip All Duplicate Files';
            end 
        %case 'Overwrite All'
            %replace everything, delete excel files
    end
end


end