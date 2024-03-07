
clc; % Clear command window.
clearvars; % Get rid of variables from prior run of this m-file.
imtool close all;  % Close all imtool figures.
format long g;
format compact;
captionFontSize = 14;
%%
baseFileName = "M116 14SEP2022 4X.tif";
%baseFileName = "M116 4X NL.tif";
%baseFileName = "MM116 4X NL.tif";
folder = fileparts(which(baseFileName)); % Determine where demo folder is (works with all versions).
fullFileName = fullfile(folder, baseFileName);
fprintf('Full File Name = "%s".\n', fullFileName);
%%
% If we get here, we should have found the image file.
originalImage = imread(fullFileName);
% Check to make sure that it is grayscale, just in case the user substituted their own image.
[rows, columns, numberOfColorChannels] = size(originalImage);
if numberOfColorChannels > 1
	promptMessage = sprintf('Your image file has %d color channels.\nThis demo was designed for grayscale images.\nDo you want me to convert it to grayscale for you so you can continue?', numberOfColorChannels);
	button = questdlg(promptMessage, 'Continue', 'Convert and Continue', 'Cancel', 'Convert and Continue');
	if strcmp(button, 'Cancel')
		fprintf(1, 'Finished running SphereImageAnalysis.m.\n');
		return;
	end
	% Do the conversion using standard book formula
	originalImage = rgb2gray(originalImage);
end

%%
% Display the grayscale image.
subplot(3, 3, 1);
imshow(originalImage);
% Maximize the figure window.
hFig1 = gcf;
hFig1.Units = 'normalized';
hFig1.WindowState = 'maximized'; % Go to full screen.
hFig1.NumberTitle = 'off'; % Get rid of "Figure 1"
hFig1.Name = 'Demo by Image Analyst'; % Put this into title bar.
% Force it to display RIGHT NOW (otherwise it might not display until it's all done, unless you've stopped at a breakpoint.)
drawnow;
caption = sprintf('Sphere Image Converted to BLK and WHT for Thresholding.');
title(caption, 'FontSize', captionFontSize);
axis('on', 'image'); % Make sure image is not artificially stretched because of screen's aspect ratio.
%%
sphere = imread(baseFileName);
sphere = rgb2gray(sphere);
sphere_imadjust = imadjust(sphere); %This seems to be the best!
sphere_histeq = histeq(sphere);
sphere_adapthisteq = adapthisteq(sphere);

montage({sphere,sphere_imadjust,sphere_histeq,sphere_adapthisteq},"Size",[1 4])
title("Original Image and Enhanced Images using imadjust, histeq, and adapthisteq")

%Choose image with highest contrast to make new original image
enhancedOriginalImage = sphere_imadjust;
subplot(3,3,2);
imshow(enhancedOriginalImage);

%%
% Just for fun, let's get its histogram and display it.
[pixelCount, grayLevels] = imhist(originalImage);
subplot(3, 3, 3);
bar(pixelCount);
title('Histogram of original image', 'FontSize', captionFontSize);
xlim([0 grayLevels(end)]); % Scale x axis manually.
grid on;

%% 


%------------------------------------------------------------------------------------------------------------------------------------------------------
% Threshold the image to get a binary image (only 0's and 1's) of class "logical."
% Method #1: using im2bw()
%   normalizedThresholdValue = 0.4; % In range 0 to 1.
%   thresholdValue = normalizedThresholdValue * max(max(originalImage)); % Gray Levels.
%   binaryImage = im2bw(originalImage, normalizedThresholdValue);       % One way to threshold to binary
% Method #2: using a logical operation.
lowThreshold = 75;
highThreshold = 175;
binaryImage = (enhancedOriginalImage >= lowThreshold) & (enhancedOriginalImage <= highThreshold);
% ========== IMPORTANT OPTION ============================================================
% Use < if you want to find dark objects instead of bright objects.
%   binaryImage = originalImage < thresholdValue; % Dark objects will be chosen if you use <.

% Do a "hole fill" to get rid of any background pixels or "holes" inside the blobs.
binaryImage = imfill(binaryImage, 'holes');

% Show the threshold as a vertical red bar on the histogram.
hold on;
maxYValue = ylim;
line([lowThreshold, lowThreshold], maxYValue, 'Color', 'r');
%line([highThreshold, highThreshold],maxYValue, 'Color','r');
% Place a text label on the bar chart showing the threshold.
annotationText1 = sprintf('Low Thresholded at %d gray levels', lowThreshold);
%annotationText2 = sprintf('High Thresholded at %d gray levels', highThreshold);
% For text(), the x and y need to be of the data class "double" so let's cast both to double.
% text(double(highThreshold + 5), double(0.5 * maxYValue(2)), annotationText1, 'FontSize', 10, 'Color', [0 .5 0]);
% text(double(highThreshold - 70), double(0.94 * maxYValue(2)), 'Background', 'FontSize', 10, 'Color', [0 0 .5]);
% text(double(highThreshold + 50), double(0.94 * maxYValue(2)), 'Foreground', 'FontSize', 10, 'Color', [0 0 .5]);

%% 
% Display the binary image.
subplot(3, 3, 4);
imshow(binaryImage);
title('Binary Image, obtained by thresholding', 'FontSize', captionFontSize);

%% Additional Processing
I  = binaryImage; 
%stdfilt
S = stdfilt(I,ones(3));
imtool(S<8); %ridge lines kinda
M = S>=8;  %map where the std is large
imtool(M);
%% another try
bw2 = ~bwareaopen(~binaryImage, 40);
subplot(3,3,5);
imshow(bw2);
%%
D = -bwdist(~binaryImage);
subplot(3,3,6);
imshow(D,[])
%%
Ld = watershed(D);
subplot(3,3,7);
imshow(label2rgb(Ld))
%%
bw2 = binaryImage;
bw2(Ld == 0) = 0;
subplot(3,3,8);
imshow(bw2)
%%
mask = imextendedmin(D,5);
subplot(3,3,9);
imshowpair(binaryImage,mask,'blend')
%%
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = binaryImage;
bw3(Ld2 == 0) = 0;
subplot(3,3,1);
imshow(bw3)
%%
props = regionprops(bw3, enhancedOriginalImage,'all');
numberOfBlobs = numel(props);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Remove connected Objects on Border
BWnobord = imclearborder(binaryImage,4);
%subplot(3,3,6);
imshow(BWnobord)
title('Cleared Border Image')
%%
bw3 = imopen(binaryImage, ones(5,5));
%subplot(3,3,5);
imshow(bw3)
title('Binary Image, adding imopen', 'FontSize', captionFontSize);
%%
bw4 = bwareaopen(bw3,40);
%subplot(3,3,6);
imshow(bw4)
title('Binary Image, adding bwareaopen', 'FontSize', captionFontSize);
%%
bw4_perim = bwperim(bw4,8);
%subplot(3,3,7);
imshow(bw4_perim)
title('Binary Image, adding bwperim', 'FontSize', captionFontSize);
%%
overlay1 = imoverlay(enhancedOriginalImage, bw4_perim, [0.3 1 0.3]);
%subplot(3,3,8);
imshow(overlay1)
title('Adjusted Binary Image, Overlayed on Original', 'FontSize', captionFontSize);
%%
mask_em = imextendedmax(enhancedOriginalImage,20);
subplot(3,3,5);
imshow(mask_em);
%%
se = strel('disk',5);
mask_em = imclose(mask_em, se); %was: ones(5,5)
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 40);
overlay2 = imoverlay(enhancedOriginalImage, bw4_perim | mask_em, [.3 1 .3]);
subplot(3,3,6);
imshow(overlay2);
%%
originalImage_c = imcomplement(enhancedOriginalImage);
subplot(3,3,7);
imshow(originalImage_c);
%%
I_mod = imimposemin(originalImage_c, ~bw4 | mask_em);
subplot(3,3,8);
imshow(I_mod);
%%
L = watershed(I_mod);
subplot(3,3,9);
imshow(label2rgb(L));
%%
props = regionprops(L, enhancedOriginalImage,'all');
numberOfBlobs = numel(props);

%% 
% Print out the measurements to the command window, and display blob numbers on the image.
textFontSize = 14;	% Used to control size of "blob number" labels put atop the image.
% Print header line in the command window.
fprintf(1,'Blob #      Mean Intensity  Area   Perimeter    Centroid       Diameter\n');
% Extract all the mean diameters into an array.
% The "diameter" is the "Equivalent Circular Diameter", which is the diameter of a circle with the same number of pixels as the blob.
% Enclosing in brackets is a nice trick to concatenate all the values from all the structure fields (every structure in the props structure array).
blobECD = [props.EquivDiameter];
% Loop over all blobs printing their measurements to the command window.
count = 0;
kStore = [];
for k = 1 : numberOfBlobs           % Loop through all blobs.
    if blobECD(k) > 100 && blobECD(k) < 2000 %minimum calculated diameter to filter out tiny detected blobs
        count = count + 1;
	    % Find the individual measurements of each blob.  They are field of each structure in the props strucutre array.
	    % You could use the bracket trick (like with blobECD above) OR you can get the value from the field of this particular structure.
	    % I'm showing you both ways and you can use the way you like best.
	    meanGL = props(k).MeanIntensity;		% Get average intensity.
	    blobArea = props(k).Area;				% Get area.
	    blobPerimeter = props(k).Perimeter;		% Get perimeter.
	    blobCentroid = props(k).Centroid;		% Get centroid one at a time
        blobDiameter(count) = blobECD(k);              %Specific diameter value
        kStore = [kStore; k];       %for indexing props after
	    % Now do the printing of this blob's measurements to the command window.
	    fprintf(1,'#%2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f\n', k, meanGL, blobArea, blobPerimeter, blobCentroid, blobDiameter(count));
	    % Put the "blob number" labels on the grayscale image that is showing the red boundaries on it.
	    text(blobCentroid(1), blobCentroid(2), num2str(count), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
newProps = props(kStore);
SphereStats.max = max(blobDiameter);
SphereStats.min = min(blobDiameter);
SphereStats.mean = mean(blobDiameter);
%%




% %%
% %TEST Prosseccing!!
% 
% 
% %use watershed analysis
% D = -bwdist(binaryImage,'euclidean');
% D(binaryImage) = -inf;  %set background to be infinitely far away
% im2 = watershed(D);
% subplot(3,3,8);
% imshow(label2rgb(im2));
% 
% 
% 
% 
% %get the region properties
% props2 = regionprops(im2,'all');
% numberOfBlobs = numel(props2);
% 
% 



%plot the results
figure;
imagesc(binaryImage);
hold on;
%plot the actual centroids
%Filter centroids
xCentroid = zeros(length(newProps),1);
yCentroid = zeros(length(newProps),1);
for ai = 1:length(newProps)
    Centroids = newProps(ai).Centroid;
    xCentroid(ai) = Centroids(1); 
    yCentroid(ai) = Centroids(2);
end
for ai = 1:length(xCentroid)
    p(1) = plot(xCentroid(ai),yCentroid(ai),'r+','markersize',10);
end
%loop over all detected objects
for ai = 1:length(newProps)
    %check if it is a circle
    if abs(newProps(ai).Perimeter-pi*newProps(ai).EquivDiameter)/newProps(ai).EquivDiameter < 0.1
        %plot the detected centroid
        p(2) = plot(newProps(ai).Centroid(1),newProps(ai).Centroid(2),'kx','markersize',10);
    end
end
legend(p,{'Actual centroid','Detected centroid'});

%%
I = imread(autofilegroup{1});
I = rgb2gray(I);
edgeIm = edge(I,"sobel",0.05);
im3 = repmat(edgeIm, [1 1 3]);
image(im3)

se90 = strel('line',1,90);
se0 = strel('line',3,0);
BWsdil = imdilate(edgeIm,[se90 se0]);
imshow(BWsdil)
title('Dilated Gradient Mask');

BWdfill = imfill(BWsdil,'holes');
imshow(BWdfill)
title('Binary Image with Filled Holes')

%Smoothing
seD = strel('diamond',1);
BWfinal = imerode(BWdfill,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal)
title('Segmented Image');
%%

%%read in your image saved from the interwebs
I = imread(autofilegroup{2});
I = rgb2gray(I);
I = imadjust(I);
S = stdfilt(I,ones(3));
imtool(S<5); %ridge lines kinda
M = S>=5;  %map where the std is large
I2 = I.*uint8(M);  %apply map
imtool(I2) %blanked out everything but ridge lines.  \
imtool(M); %map
%%
%Dilate
se90 = strel('line',5,90);
se0 = strel('line',5,0);
BWsdil = imdilate(M,[se90 se0]);
imshow(BWsdil)
title('Dilated Gradient Mask');
%%
bw2 = bwareaopen(BWsdil, 40);
imshow(bw2);
%%
bw2 = imfill(bw2,'holes');
imshow(bw2);
title('Binary Image with Filled Holes')
%%
%Smoothing
seD = strel('diamond',3);
BWfinal = imerode(bw2,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal)
title('Segmented Image');
%%
D = -bwdist(~BWfinal);
imtool(D,[])
%%
Ld = watershed(D);
imtool(label2rgb(Ld))
%%
bw2 = BWfinal;
bw2(Ld == 0) = 0;
imtool(bw2)
%%
mask = imextendedmin(D,1);
imshowpair(BWfinal,mask,'blend')
%%
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = BWsdil;
bw3(Ld2 == 0) = 0;
imshow(bw3)
%%
props = regionprops(bw3, I,'all');
numberOfBlobs = numel(props);




%%
I = imread(autofilegroup{2});
I = rgb2gray(I);
I = imadjust(I);
imshow(I);
%% Another Test
r = drawrectangle;
%%
%mask = createMask(r);
bw = activecontour(T,mask,1000,'edge');
hold on;
visboundaries( bw,'Color','r'); 
%%
[~,threshold] = edge(I,'sobel');
fudgeFactor = 0.5;
BWs = edge(I,'sobel',threshold * fudgeFactor);
imshow(BWs)
%%
S = stdfilt(bw2,ones(3));
imshow(S<5); %ridge lines kinda
M = S>=5;  %map where the std is large
I2 = bw2.*uint8(M);  %apply map
imtool(I2) %blanked out everything but ridge lines.  \
imtool(M); %map
%%
%Dilate
se90 = strel('line',10,90);
se0 = strel('line',10,0);
BWsdil = imdilate(BWs,[se90 se0]);
imshow(BWsdil)
title('Dilated Gradient Mask');

%%
BWfill = imfill(BWsdil,'holes');
imshow(BWfill);
title('Binary Image with Filled Holes')

%%
seD = strel('diamond',1);
BWfinal = imerode(BWfill,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal)
title('Segmented Image');


% function [BW,maskedImage] = detectCircles(RGB)
% %segmentImage Segment image using auto-generated code from imageSegmenter app
% %  [BW,MASKEDIMAGE] = segmentImage(RGB) segments image RGB using
% %  auto-generated code from the imageSegmenter app. The final segmentation
% %  is returned in BW, and a masked image is returned in MASKEDIMAGE.
% 
% % Auto-generated by imageSegmenter app on 20-Jul-2022
% %----------------------------------------------------
% 
% 
% % Convert RGB image into L*a*b* color space.
% X = rgb2lab(RGB);
% 
% % Find circles
% [centers,radii,~] = imfindcircles(RGB,[100 500],'ObjectPolarity','dark','Sensitivity',0.98);
% BW = false(size(RGB,1),size(RGB,2));
% [Xgrid,Ygrid] = meshgrid(1:size(BW,2),1:size(BW,1));
% for n = 1:15
%     BW = BW | (hypot(Xgrid-centers(n,1),Ygrid-centers(n,2)) <= radii(n));
% end
% 
% % Create masked image.
% maskedImage = RGB;
% maskedImage(repmat(~BW,[1 1 3])) = 0;
% end
%%
%Attempt 30SEP2022
I = imread(autofilegroup{1});
I = rgb2gray(I);
I = imadjust(I);

B_thresh=graythresh(I);
BW = imbinarize(I,B_thresh);
b1 = imfill(im2bw(I,0.9), 'holes');
% BW = imbinarize(B_imadjust);
% BW = imbinarize(B_imadjust,B_thresh);
imshow(BW)
%%
B_w = edge(BW,'canny',0.86);
figure,imshow(B_w)
title('Edge detection')
%%
% Simple segmentation
tic
mask = zeros(size(I));
mask(25:end-25,25:end-25) = 1;
bw = activecontour(I,mask,20000);
figure;
imshow(bw)
toc
%%
BW2 = bwareaopen(bw, 700);
figure,imshow(BW2)
title('segmentation')
% 
% 
% 
% 
% 
% 

%%
FullFileName = "C:\Users\StevenSummey\Documents\MATLAB\Image Processing\Microscope Images\M116 14SEP2022 4X.tif";
im = imread(FullFileName);
load('C:\Users\StevenSummey\Documents\MATLAB\Image Processing\Microscope Images\hNP12SEP2022B-1 M-PROCESSED\M116-14SEP2022.mat')
im = rgb2gray(im);
e = edge(im, 'canny');
imshow(e);
%%
FullFileName = "C:\Users\StevenSummey\Documents\MATLAB\Image Processing\Microscope Images\M116 14SEP2022 4X.tif";
im = imread(FullFileName);
im = rgb2gray(im);
BW = imbinarize(im);
BW = edge(BW, 'canny');
[B,L,n] = bwboundaries(BW,'noholes');
imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
%%
im = imread(FullFileName);
im = rgb2gray(im);
BW = imbinarize(im);
% BW = bwareaopen(BW, 40);
BW = edge(BW, 'sobel');
% BW2 = bwmorph(BW,'remove',Inf);
imshow(BW)
[B,L,n] = bwboundaries(BW,'noholes');
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1)
end
%%
tic
[accum, circen, cirrad] = CircularHough_Grd(im, [25 500]);
toc
%%
figure(1); imagesc(accum); axis image;
 title('Accumulation Array from Circular Hough Transform');
 figure(2); imagesc(im); colormap('gray'); axis image;
 hold on;
 plot(circen(:,1), circen(:,2), 'r+');
 for k = 1 : size(circen, 1)
     DrawCircle(circen(k,1), circen(k,2), cirrad(k), 32, 'b-');
 end
 hold off;
 title(['Raw Image with Circles Detected ', ...
     '(center positions and radii marked)']);
 figure(3); surf(accum, 'EdgeColor', 'none'); axis ij;
 title('3-D View of the Accumulation Array');
%%
%%read in your image saved from the interwebs
I = imread(baseFileName);
I = rgb2gray(I);
I = imadjust(I);
S = stdfilt(I,ones(3));
imtool(S<10); %ridge lines kinda
M = S>=10;  %map where the std is large
I2 = I.*uint8(M);  %apply map
imtool(I2) %blanked out everything but ridge lines.  \
imtool(M); %map
%%
BW = imclearborder(M, 8);
imshow(BW);
BW2 = bwareaopen(BW, 200); 
imshow(BW2);
%%
BWsdil = BW2;
dilateState = 1;
thickness = 0;
while dilateState == 1
    button = menu('Target Sphere Edges Connected?', 'YES','NO');
    switch button
        case 1
            dilateState = 0; 
        case 2
            thickness = thickness + 5; 
            se90 = strel('line',thickness,90);
            se0 = strel('line',thickness,0);
            BWsdil = imdilate(BW2,[se90 se0]);
            imshow(BWsdil);
            title('Dilated Gradient Mask');
            dilateState = 1; 
    end
end
%%
%Dilate
se90 = strel('line',20,90);
se0 = strel('line',20,0);
BWsdil = imdilate(BW2,[se90 se0]);
imshow(BWsdil)
title('Dilated Gradient Mask');

%%
BWfill = imfill(BWsdil,'holes');
imshow(BWfill);
title('Binary Image with Filled Holes')
%%
seD = strel('diamond',5);
BWfinal = imerode(BWfill,seD);
BWfinal = imerode(BWfinal,seD);
imshow(BWfinal)
title('Segmented Image');
%%
% Simple segmentation
tic
mask = zeros(size(I2));
mask(25:end-25,25:end-25) = 1;
bw = activecontour(I2,mask,5000);
figure;
imshow(bw)
toc
%%
A = imread(baseFileName);
A = rgb2gray(A);
B = imbinarize(A,'global');
se = strel('disk',5);
C = imclose(~B,se);
imshow(C)
%%
%%read in your image saved from the interwebs
I = imread(baseFileName);
I = rgb2gray(I);
I = imadjust(I);
S = stdfilt(I,ones(3));
imshow(S<5); %ridge lines kinda
M = S>=5;  %map where the std is large
I2 = I.*uint8(M);  %apply map
imshow(I2) %blanked out everything but ridge lines.  \
imshow(M); %map
%%
sM = bwmorph(M,'hbreak',Inf);
imshow(sM)
%%
se = strel('disk',2);
C = bwmorph(sM,'close',Inf);
%C = imclose(sM,se);
imshow(C)
%%
D = imfill(C, 'holes');
imshow(D)
%%
D = -bwdist(~C);
imtool(D,[])
%%
Ld = watershed(D);
imtool(label2rgb(Ld))
%%
bw2 = C;
bw2(Ld == 0) = 0;
imtool(bw2)
%%
mask = imextendedmin(D,2);
imshowpair(C,mask,'blend')
%%
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = C;
bw3(Ld2 == 0) = 0;
imshow(bw3)
%%
props = regionprops(bw3, I,'all');
numberOfBlobs = numel(props);
%%

%%read in your image saved from the interwebs
I = imread(baseFileName);
I = rgb2gray(I);
I = imadjust(I);
S = stdfilt(I,ones(3));
imshow(S<10); %ridge lines kinda
M = S>=10;  %map where the std is large
% M = ~M;
I2 = I.*uint8(M);  %apply map
imshow(I2) %blanked out everything but ridge lines.  \
imshow(M);
%%
bw1 = edge(I2,'canny');
imshow(bw1)
%%
bw1 = bwareaopen(bw1, 200);
imshow(bw1)
%%

%%
bw2 = imfill(bw1,'holes');
imshow(bw2)






