function [SphereStats,autoStore,percentDiff,newProps] = imagePropertyAnalysis(finalImage,originalImage,autoStore,...
    referenceStats, SphereStats, distancePerPixel,init)

    if init == 0
        props = regionprops(finalImage,originalImage, 'all');
        numberOfSpheres = numel(props);
        % Print out the measurements to the command window, and display blob numbers on the image.
        % Print header line in the command window.
        fprintf(1,'Sphere #      Mean Intensity  Area   Perimeter    Centroid       Diameter    Circularity\n');
        % Extract all the mean diameters into an array.
        % The "diameter" is the "Equivalent Circular Diameter", which is the diameter of a circle with the same number of pixels as the blob.
        % Enclosing in brackets is a nice trick to concatenate all the values from all the structure fields (every structure in the props structure array).
        blobECD = [props.EquivDiameter];
        blobCirc = [props.Circularity];
        circThreshold = 0.80; %more likely to include oblong spheres
        minThreshold = referenceStats.min/distancePerPixel; %minimum pixel threshold, constant per batch
        maxThreshold = max(blobECD)*1.1; %will change with each individual analysis, scale up to catch higher value
    
        %%% AUTO STORE VALUE %%%
        autoStore = [autoStore string(minThreshold)];
        autoStore = [autoStore string(maxThreshold)];
        autoStore = [autoStore "AUTO"]; %final tag
        autoStore = strjoin(autoStore,'_');
        %%% ---------------- %%%
        % Loop over all blobs printing their measurements to the command window.
        count = 0;
        kStore = [];
        for k = 1 : numberOfSpheres          % Loop through all blobs.
            if blobECD(k) >= minThreshold && blobECD(k) <= maxThreshold && blobCirc(k) >= circThreshold%minimum calculated diameter to filter out tiny detected blobs
                count = count + 1;
                % Find the individual measurements of each blob.  They are field of each structure in the props strucutre array.
                % You could use the bracket trick (like with blobECD above) OR you can get the value from the field of this particular structure.
                % I'm showing you both ways and you can use the way you like best.
                meanGL = props(k).MeanIntensity;		% Get average intensity.
                blobArea = props(k).Area;				% Get area.
                blobPerimeter = props(k).Perimeter;		% Get perimeter.
                blobCentroid = props(k).Centroid;		% Get centroid one at a time
                blobDiameter(count) = blobECD(k);              %Specific diameter value
                blobCircularity(count) = blobCirc(k); 
                kStore = [kStore; k];       %for indexing props after
                % Now do the printing of this blob's measurements to the command window.
                fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f %8.1f %8.2f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobDiameter(count), blobCircularity(count));
                % Put the "blob number" labels on the grayscale image that is showing the red boundaries on it.
                %text(blobCentroid(1), blobCentroid(2), num2str(count), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            end
        end
        %Convert applicable fields to microns using calibration value
        newProps = props(kStore);
        sphereDiameter = zeros(numel(newProps),1);
        circularity = zeros(numel(newProps),1);
        solidity = zeros(numel(newProps),1);
        eccentricity = zeros(numel(newProps),1);
        for k = 1:numel(newProps)
            newProps(k).EquivDiameter = newProps(k).EquivDiameter.*distancePerPixel;
            sphereDiameter(k) = newProps(k).EquivDiameter; 
            circularity(k) = newProps(k).Circularity;
            solidity(k) = newProps(k).Solidity; 
            eccentricity(k) = newProps(k).Eccentricity; 
        end
        numberOfSpheres = numel(newProps);
        sprintf('%d spheres ultimately detected.\n', numberOfSpheres);
        SphereStats.count = numberOfSpheres;
        SphereStats.max = max(sphereDiameter); % converted
        SphereStats.min = min(sphereDiameter); % converted
        SphereStats.mean = mean(sphereDiameter); % converted
    
        percentDiff = (abs(referenceStats.mean - SphereStats.mean)/referenceStats.mean)*100;
    
        SphereStats.median = median(sphereDiameter); %converted
        SphereStats.std = std(sphereDiameter); %converted
        %find mode bin size
        [N, edges] = histcounts(sphereDiameter,20);
        maxIndex = find(N==max(N));
        modeBin = [];
        if numel(maxIndex) > 1
            edgeStore = zeros(numel(maxIndex)*2, 1);
            for j = 1:numel(maxIndex)
                edgeStore(j) = edges(maxIndex(j));
                edgeStore(j+1) = edges(maxIndex(j)+1);
                if j == 1                         
                    modeBin = [modeBin [num2str(edgeStore(j)) '-' num2str(edgeStore(j+1))]];
                else
                    modeBin = [modeBin ', ' [num2str(edgeStore(j)) '-' num2str(edgeStore(j+1))]];
                end
            end
            SphereStats.modeBin = modeBin;
        else
            edge1 = edges(maxIndex);
            edge2 = edges(maxIndex+1);
            SphereStats.modeBin = [num2str(edge1) '-' num2str(edge2)];
        end
        SphereStats.AverageCircularity = mean(circularity);
        SphereStats.AverageSolidity = mean(solidity);
        SphereStats.AverageEccentricity = mean(eccentricity);
        SphereStats.PolyDispersity = (SphereStats.std/SphereStats.mean)^2; 
        SphereStats.FLAG = 0;
        SphereStats.Comments = "Batch Processing";

    elseif init == 1 %DLNN settting --> incorporate circularity filter

        props = regionprops(finalImage,originalImage, 'all');
        numberOfSpheres = numel(props);
        % Print out the measurements to the command window, and display blob numbers on the image.
        % Print header line in the command window.
        fprintf(1,'Sphere #      Mean Intensity  Area   Perimeter    Centroid       Diameter    Circularity\n');
        % Extract all the mean diameters into an array.
        % The "diameter" is the "Equivalent Circular Diameter", which is the diameter of a circle with the same number of pixels as the blob.
        % Enclosing in brackets is a nice trick to concatenate all the values from all the structure fields (every structure in the props structure array).
        blobECD = [props.EquivDiameter];
        blobCirc = [props.Circularity];
        circThreshold = 0.80; %more likely to include oblong spheres
        
        %rmoutliers with bias toward bottom threshold (excess segmentation
        %buffer)
        [~,~,~,L,~,~] = rmoutliers(blobECD,"percentiles",[25 100]);
        % for later use [B,TFrm,TFoutlier,L,U,C]
        % plot(A)
        % hold on
        % plot(find(~TFrm),B,"o-")
        % yline([L U C],":",["Lower Threshold","Upper Threshold","Center Value"])
        % legend("Original Data","Cleaned Data")
        %minThreshold = mean(blobECD)-2*std(blobECD); %minimum pixel threshold, three standard deviations from mean
        minThreshold = L; 
        if minThreshold < 30
            minThreshold = 30; 
        end
        maxThreshold = max(blobECD)*1.1; %will change with each individual analysis
        if minThreshold == 30 && maxThreshold > 200 %run through again to avoid debris
            count = 0;
            kStore = [];
            for k = 1 : numberOfSpheres          % Loop through all blobs.
                %size and circularity filter
                if blobECD(k) >= minThreshold && blobECD(k) <= maxThreshold && blobCirc(k) >= circThreshold
                    %minimum calculated diameter to filter out tiny detected blobs
                    count = count + 1;
                    % Find the individual measurements of each blob.  They are field of each structure in the props strucutre array.
                    % You could use the bracket trick (like with blobECD above) OR you can get the value from the field of this particular structure.
                    % I'm showing you both ways and you can use the way you like best.
                    meanGL = props(k).MeanIntensity;		% Get average intensity.
                    blobArea = props(k).Area;				% Get area.
                    blobPerimeter = props(k).Perimeter;		% Get perimeter.
                    blobCentroid = props(k).Centroid;		% Get centroid one at a time
                    blobDiameter(count) = blobECD(k);              %Specific diameter value
                    blobCircularity(count) = blobCirc(k); 
                    kStore = [kStore; k];       %for indexing props after
                    % Now do the printing of this blob's measurements to the command window.
                    fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f %8.1f %8.2f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobDiameter(count), blobCircularity(count));
                    % Put the "blob number" labels on the grayscale image that is showing the red boundaries on it.
                    %text(blobCentroid(1), blobCentroid(2), num2str(count), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end

            props = props(kStore);
            numberOfSpheres = numel(props);
            fprintf(1,'Sphere #      Mean Intensity  Area   Perimeter    Centroid       Diameter    Circularity\n');
            blobECD = [props.EquivDiameter];
            blobCirc = [props.Circularity];
            circThreshold = 0.80; %more likely to include oblong spheres
            
            %rmoutliers with bias toward bottom threshold (excess segmentation
            %buffer)
            [~,~,~,L,~,~] = rmoutliers(blobECD,"percentiles",[25 100]);
            minThreshold = L; 
            if minThreshold < 50
                minThreshold = 50; 
            end
            maxThreshold = max(blobECD)*1.1; %will change with each individual analysis
            %%% AUTO STORE VALUE %%%
            autoStore = [autoStore string(minThreshold)];
            autoStore = [autoStore string(maxThreshold)];
            autoStore = [autoStore "AUTO"]; %final tag
            autoStore = strjoin(autoStore,'_');
            %%% ---------------- %%%
            count = 0;
            kStore = [];
            for k = 1 : numberOfSpheres          % Loop through all blobs.
                %size and circularity filter
                if blobECD(k) >= minThreshold && blobECD(k) <= maxThreshold && blobCirc(k) >= circThreshold
                    %minimum calculated diameter to filter out tiny detected blobs
                    count = count + 1;
                    % Find the individual measurements of each blob.  They are field of each structure in the props strucutre array.
                    % You could use the bracket trick (like with blobECD above) OR you can get the value from the field of this particular structure.
                    % I'm showing you both ways and you can use the way you like best.
                    meanGL = props(k).MeanIntensity;		% Get average intensity.
                    blobArea = props(k).Area;				% Get area.
                    blobPerimeter = props(k).Perimeter;		% Get perimeter.
                    blobCentroid = props(k).Centroid;		% Get centroid one at a time
                    blobDiameter(count) = blobECD(k);              %Specific diameter value
                    blobCircularity(count) = blobCirc(k); 
                    kStore = [kStore; k];       %for indexing props after
                    % Now do the printing of this blob's measurements to the command window.
                    fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f %8.1f %8.2f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobDiameter(count), blobCircularity(count));
                    % Put the "blob number" labels on the grayscale image that is showing the red boundaries on it.
                    %text(blobCentroid(1), blobCentroid(2), num2str(count), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end
        else
            %%% AUTO STORE VALUE %%%
            autoStore = [autoStore string(minThreshold)];
            autoStore = [autoStore string(maxThreshold)];
            autoStore = [autoStore "AUTO"]; %final tag
            autoStore = strjoin(autoStore,'_');
            %%% ---------------- %%%
            % Loop over all blobs printing their measurements to the command window.
            count = 0;
            kStore = [];
            for k = 1 : numberOfSpheres          % Loop through all blobs.
                %size and circularity filter
                if blobECD(k) >= minThreshold && blobECD(k) <= maxThreshold && blobCirc(k) >= circThreshold
                    %minimum calculated diameter to filter out tiny detected blobs
                    count = count + 1;
                    % Find the individual measurements of each blob.  They are field of each structure in the props strucutre array.
                    % You could use the bracket trick (like with blobECD above) OR you can get the value from the field of this particular structure.
                    % I'm showing you both ways and you can use the way you like best.
                    meanGL = props(k).MeanIntensity;		% Get average intensity.
                    blobArea = props(k).Area;				% Get area.
                    blobPerimeter = props(k).Perimeter;		% Get perimeter.
                    blobCentroid = props(k).Centroid;		% Get centroid one at a time
                    blobDiameter(count) = blobECD(k);              %Specific diameter value
                    blobCircularity(count) = blobCirc(k); 
                    kStore = [kStore; k];       %for indexing props after
                    % Now do the printing of this blob's measurements to the command window.
                    fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f %8.1f %8.2f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobDiameter(count), blobCircularity(count));
                    % Put the "blob number" labels on the grayscale image that is showing the red boundaries on it.
                    %text(blobCentroid(1), blobCentroid(2), num2str(count), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
                end
            end
        end
        %Convert applicable fields to microns using calibration value
        newProps = props(kStore);
        assignin('base',"newProps",newProps);
        sphereDiameter = zeros(numel(newProps),1);
        circularity = zeros(numel(newProps),1);
        solidity = zeros(numel(newProps),1);
        eccentricity = zeros(numel(newProps),1);
        for k = 1:numel(newProps)
            newProps(k).EquivDiameter = newProps(k).EquivDiameter.*distancePerPixel;
            sphereDiameter(k) = newProps(k).EquivDiameter; 
            circularity(k) = newProps(k).Circularity;
            solidity(k) = newProps(k).Solidity; 
            eccentricity(k) = newProps(k).Eccentricity; 
        end
        numberOfSpheres = numel(newProps);
        sprintf('%d spheres ultimately detected.\n', numberOfSpheres);
        SphereStats.count = numberOfSpheres;
        SphereStats.max = max(sphereDiameter); % converted
        SphereStats.min = min(sphereDiameter); % converted
        SphereStats.mean = mean(sphereDiameter); % converted
           
        SphereStats.median = median(sphereDiameter); %converted
        SphereStats.std = std(sphereDiameter); %converted
        SphereStats.percCV = (SphereStats.std/SphereStats.mean)*100; 
        %find mode bin size
        [N, edges] = histcounts(sphereDiameter,20);
        maxIndex = find(N==max(N));
        modeBin = [];
        if numel(maxIndex) > 1
            edgeStore = zeros(numel(maxIndex)*2, 1);
            for j = 1:numel(maxIndex)
                edgeStore(j) = edges(maxIndex(j));
                edgeStore(j+1) = edges(maxIndex(j)+1);
                if j == 1                         
                    modeBin = [modeBin [num2str(edgeStore(j)) '-' num2str(edgeStore(j+1))]];
                else
                    modeBin = [modeBin ', ' [num2str(edgeStore(j)) '-' num2str(edgeStore(j+1))]];
                end
            end
            SphereStats.modeBin = modeBin;
        else
            edge1 = edges(maxIndex);
            edge2 = edges(maxIndex+1);
            SphereStats.modeBin = [num2str(edge1) '-' num2str(edge2)];
        end
        SphereStats.AverageCircularity = mean(circularity);
        SphereStats.AverageSolidity = mean(solidity);
        SphereStats.AverageEccentricity = mean(eccentricity);
        SphereStats.PolyDispersity = (SphereStats.std/SphereStats.mean)^2; 
        if SphereStats.count <= 1
            SphereStats.FLAG = 1; 
            SphereStats.Comments = "Not enough sphere counts";
        else
            SphereStats.FLAG = 0;
            SphereStats.Comments = "Batch Processing-DLNN";
        end
        

        percentDiff = 0; %not used in this init state
    end
end