%defines the smallest image size the network can process
inputSize = [32 32 3];
imgLayer = imageInputLayer(inputSize);

%% Create Downsampling Network
% convolution and ReLu layershe convolution layer padding is selected such 
% %that the output size of the convolution layer is the same as the input size. 
% This makes it easier to construct a network because the input and output sizes 
% between most layers remain the same as you progress through the network.
filterSize = 3;
numFilters = 32;
conv = convolution2dLayer(filterSize,numFilters,'Padding',1);
relu = reluLayer();

% The downsampling is performed using a max pooling layer. Create a max 
% pooling layer to downsample the input by a factor of 2 by setting the 
% 'Stride' parameter to 2.
poolSize = 2;
maxPoolDownsample2x = maxPooling2dLayer(poolSize,'Stride',2);

%Stack the convolution, ReLU, and max pooling layers to create a network 
% that downsamples its input by a factor of 4.
downsamplingLayers = [
    conv
    relu
    maxPoolDownsample2x
    conv
    relu
    maxPoolDownsample2x
    ];

%% Create Upsampling Network
%The upsampling is done using the tranposed convolution layer (also commonly 
% referred to as "deconv" or "deconvolution" layer). When a transposed 
% convolution is used for upsampling, it performs the upsampling and the 
% filtering at the same time.

%Create a transposed convolution layer to upsample by 2. 
filterSize = 4;
transposedConvUpsample2x = transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping',1);
%The 'Cropping' parameter is set to 1 to make the output size equal twice the input size.

%Stack the transposed convolution and relu layers. An input to this set of layers is upsampled by 4.
upsamplingLayers = [
    transposedConvUpsample2x
    relu
    transposedConvUpsample2x
    relu
    ];

%% Create a Pixel Classification Layer
%The final set of layers are responsible for making pixel classifications. 
% These final layers process an input that has the same spatial dimensions 
% (height and width) as the input image. However, the number of channels (third dimension) 
% is larger and is equal to number of filters in the last transposed convolution layer. 
% This third dimension needs to be squeezed down to the number of classes we wish to segment. 
% This can be done using a 1-by-1 convolution layer whose number of filters equal the number of classes, e.g. 3.

%Create a convolution layer to combine the third dimension of the input 
% feature maps down to the number of classes.
numClasses = 3;
conv1x1 = convolution2dLayer(1,numClasses);

%Following this 1-by-1 convolution layer are the softmax and pixel classification layers. 
% These two layers combine to predict the categorical label for each image pixel.
finalLayers = [
    conv1x1
    softmaxLayer()
    pixelClassificationLayer()
    ];

%% Stack all Layers
%Stack all the layers to complete the semantic segmentation network. 
net = [
    imgLayer    
    downsamplingLayers
    upsamplingLayers
    finalLayers
    ];
%NETWORK IS READY TO BE TRAINED USING trainNetwork from DLT

%% Train Semantic Segmentation Network
% %Load the training data
% dataSetDir = fullfile(pwd,'sphereTrainingData2/GroundTruthProject');
% imageDir = fullfile(dataSetDir, 'signalInfo.mat');
% S = load(imageDir);
% imagePaths = S.S.SignalSources; 
% %imageLogDir = fullfile(dataSetDir, 'pixelLabelFlags.mat');
% %R = load(imageLogDir);
% %imageLogicals = R.S; 
% %imageDir_new = imagePaths(imageLogicals);
% labelDir = fullfile(dataSetDir,'PixelLabelData');

% when gTruth is exported

%Create an image datastore for the images.
imds = gTruth.DataSource.Source;
pxDir = gTruth.LabelData.PixelLabelData;

%Create a pixelLabelDatastore for the ground truth pixel labels.
classNames = ["Sphere", "Debris", "Background"];
labelIDs = [1 2 3];
pxds = pixelLabelDatastore(pxDir,classNames, labelIDs);

%Visualize training images and ground truth pixel labels.
I = read(imds);
C = read(pxds);

I = imresize(I,5);
L = imresize(uint8(C{1}),5);
imshowpair(I,L,'montage')


%%
%Create a semantic segmentation network. This network uses a simple 
% semantic segmentation network based on a downsampling and upsampling design. 

%resize images and pixels
targetSize = [1944 2592 3]; %smallest image dimension
% auimds = augmentedImageDatastore(targetSize, imds);
% pxdsReSz = transform(pxds,@(x) imresize(x,targetSize));
%%
numFilters = 64;
filterSize = 3;
numClasses = 3;
layers = [
    imageInputLayer(targetSize)
    convolution2dLayer(filterSize,numFilters,'Padding',1)
    reluLayer()
    maxPooling2dLayer(2,'Stride',2)
    convolution2dLayer(filterSize,numFilters,'Padding',1)
    reluLayer()
    transposedConv2dLayer(4,numFilters,'Stride',2,'Cropping',1);
    convolution2dLayer(1,numClasses);
    softmaxLayer()
    pixelClassificationLayer()
    ];

%Setup training options.
opts = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod', 5,...
    'LearnRateDropFactor', 0.2,...
    'InitialLearnRate',0.01, ...
    'Momentum',0.9,...
    'MaxEpochs',50, ...
    'MiniBatchSize', 4,...
    'Shuffle','every-epoch',...
    'Verbose',1,...
    'VerboseFrequency',1,...
    'CheckpointPath','C:\Users\StevenSummey\Documents\MATLAB\Image Processing\DeepLearningCheckpoint',...
    'CheckpointFrequency',5,...
    'CheckpointFrequencyUnit','epoch',...
    'Plots','training-progress');
 
%Combine the image and pixel label datastore for training
trainingData = combine(imds,pxds);

%Train the network.
net = trainNetwork(trainingData,layers,opts);

%%
%gTruth saved in sphereTrainingData3
gTruth = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData3\gTruth.mat');
imds_1 = gTruth.gTruth.DataSource.Source.Files;
%gTruth2 needs to have image paths changed
gTruth2 = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData4\gTruth2.mat');
pathPre = 'C:\Users\StevenSummey\Documents\MATLAB\Microscope Images\Testing';
pathPost = [pwd '\Microscope Images\Training'];
alternativePaths = [string(pathPre), string(pathPost)];
unresolvedPaths = changeFilePaths(gTruth2.gTruth, alternativePaths);
alternativePaths = [string(pathPost), string(pathPre)];

%gTruth3 needs to have image paths changed
gTruth3 = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData5\gTruth3.mat');
pathPre = 'C:\Users\StevenSummey\Documents\MATLAB\Microscope Images\TESTING_LIGHT';
pathPost = [pwd '\Microscope Images\Training'];
alternativePaths = [string(pathPre), string(pathPost)];
unresolvedPaths = changeFilePaths(gTruth3.gTruth3, alternativePaths);
alternativePaths = [string(pathPost), string(pathPre)];

%gTruth4 needs to have image paths changed
gTruth4 = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData7\gTruth4.mat');
pathPre = 'C:\Users\StevenSummey\Documents\MATLAB\Microscope Images\Testing_NO LIGHT';
pathPost = [pwd '\Microscope Images\Training'];
alternativePaths = [string(pathPre), string(pathPost)];
unresolvedPaths = changeFilePaths(gTruth4.gTruth, alternativePaths);
alternativePaths = [string(pathPost), string(pathPre)];

imds_2 = gTruth2.gTruth.DataSource.Source.Files;
imds_3 = gTruth3.gTruth3.DataSource.Source.Files; 
imds_4 = gTruth4.gTruth.DataSource.Source.Files; 
imds_combined = [imds_1; imds_2; imds_3; imds_4];
imds = imageDatastore(imds_combined);
imds.ReadFcn = @customReadDatastoreImage;

pxDir = gTruth.gTruth.LabelData.PixelLabelData;
pxDir2 = gTruth2.gTruth.LabelData.PixelLabelData;
pxDir3 = gTruth3.gTruth3.LabelData.PixelLabelData;
pxDir4 = gTruth4.gTruth.LabelData.PixelLabelData; 
pxDir_combined = [pxDir; pxDir2; pxDir3; pxDir4];

classNames = ["Sphere", "Debris", "Background"];
labelIDs = [1 2 3];
pxds = pixelLabelDatastore(pxDir_combined,classNames, labelIDs);
pxds.ReadFcn = @customReadDatastoreLabel;

%targetSize = [1944 2592 3]; %smallest image dimension [972 1296 3], [486 648 3] [1944 2592 3], [1200 1600 3]


%% PARSE LIGHT AND NO LIGHT IMAGES
%Run script above for aggregated files and labels
%TrainingLIGHT is imported data store
load('TrainingLIGHT.mat')
TrainingLIGHT.ReadFcn = @customReadDatastoreImage;
A = cellfun(@(x) strsplit(x,'\'),TrainingLIGHT.Files,'UniformOutput',false);
B = cellfun(@(x) x{end},A,'UniformOutput',false);
C = cellfun(@(x) strrep(x,'.tif','.png'), B,'UniformOutput',false);
%D = cellfun(@(x) x{1:end-1},C,'UniformOutput',false);
% TrainingLIGHT_pxds = strings(length(TrainingLIGHT.Files),1);
pxDir_search = pxDir_combined; 
for i = 1:length(C)
    for j = 1:length(pxDir_search)
        if contains(string(pxDir_search{j}), C(i))
            TrainingLIGHT_pxds(i,1) = pxDir_search(j);
            pxDir_search(j) = {'~'};
            break; 
        end
    end
end

load('TrainingNOLIGHT.mat')
TrainingNOLIGHT.ReadFcn = @customReadDatastoreImage;
A = cellfun(@(x) strsplit(x,'\'),TrainingNOLIGHT.Files,'UniformOutput',false);
B = cellfun(@(x) x{end},A,'UniformOutput',false);
C = cellfun(@(x) strrep(x,'.tif','.png'), B,'UniformOutput',false);
%D = cellfun(@(x) x{1:end-1},C,'UniformOutput',false);
% TrainingLIGHT_pxds = strings(length(TrainingLIGHT.Files),1);
pxDir_search = pxDir_combined; 
for i = 1:length(C)
    for j = 1:length(pxDir_search)
        if contains(string(pxDir_search{j}), C(i))
            TrainingNOLIGHT_pxds(i,1) = pxDir_search(j);
            pxDir_search(j) = {'~'};
            break; 
        end
    end
end
%% Organize data again
imds_LIGHT = TrainingLIGHT; %already has custom function
imds_NOLIGHT = TrainingNOLIGHT; 

classNames = ["Sphere", "Debris", "Background"];
labelIDs = [1 2 3];
pxds_LIGHT = pixelLabelDatastore(TrainingLIGHT_pxds,classNames, labelIDs);
pxds_LIGHT.ReadFcn = @customReadDatastoreLabel;%%

pxds_NOLIGHT = pixelLabelDatastore(TrainingNOLIGHT_pxds,classNames, labelIDs);
pxds_NOLIGHT.ReadFcn = @customReadDatastoreLabel;%%

%% CREATE NEW GROUNDTRUTH LABEL
%resize images and labels
targetSize = [1944 2592];
for i = 1:numel(TrainingNOLIGHT.Files)
    im_file = strsplit(TrainingNOLIGHT.Files{i},'\');
    im_file = im_file{end};
    px_file = strsplit(pxds_NOLIGHT.Files{i},'\');
    px_file = px_file{end};
    I = readimage(TrainingNOLIGHT,i);
    L = readimage(pxds_NOLIGHT,i); % added lines: 
    classes = ["Sphere","Background"]; %adjusted for no debris "Debris" , 
    ids = [1,3]; %adjusted for no debris ,2
    %C = categorical(L,ids,classes);
    I_new = imresize(I,targetSize);
    L_new = imresize(L,targetSize); % [972 1296], [486 648], [1944 2592]
    L_new(L_new == "Debris") = "Background"; %change debris to background
    cl = categorical(classes);
    [~,idx] = ismember(L_new, cl);
    search = find(ismember(idx, 0) == 1);
    if isempty(search)
        L_new = uint8(ids(idx));
    else
        idx(search) = 3; %background
        L_new = uint8(ids(idx));
    end
    imwrite(I_new,['Microscope Images\TrainingNOLIGHT_NoDebris_1944x2592\Images\' num2str(i) '-' im_file]);
    imwrite(L_new,['Microscope Images\TrainingNOLIGHT_NoDebris_1944x2592\PixelLabels\'  num2str(i) '-' px_file]);
end
%%
imDir = fullfile([pwd '\Microscope Images\TrainingNOLIGHT_NoDebris_1944x2592\Images']);
imDirFiles = dir(fullfile(imDir,'*.tif'));
for i = 1:numel(imDirFiles)
    fullFile{i,1} = [imDirFiles(i).folder '\' imDirFiles(i).name];
end

imdsNOLIGHT = imageDatastore(fullFile); 
%%
pxDir = fullfile([pwd '\Microscope Images\TrainingNOLIGHT_NoDebris_1944x2592\PixelLabels']);
pxDirFiles = dir(fullfile(pxDir,'*.png'));
for i = 1:numel(pxDirFiles)
    fullFile{i,1} = [pxDirFiles(i).folder '\' pxDirFiles(i).name];
end
classNames = ["Sphere", "Background"];
labelIDs = [1 3];
pxdsNOLIGHT = pixelLabelDatastore(fullFile,classNames,labelIDs); 
%%
load('resizedImages_NL_NoDebris_1944x2592.mat') %imdsLIGHT
imdsNOLIGHT.ReadFcn = @customReadDatastoreImage;
load('resizedPixelLabels_NL_NoDebris_1944x2592.mat') %pxdsLIGHT
pxdsNOLIGHT.ReadFcn = @customReadDatastoreLabel;
dataSource = groundTruthDataSource(imdsNOLIGHT.Files);
load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData6\GroundTruthProject\labelDefinitions.mat') %S
labelData = table(pxdsNOLIGHT.Files,'VariableNames',{'PixelLabelData'});
gTruth = groundTruth(dataSource,S,labelData);

imageLabeler(gTruth)
%% No Debris TEST
%gTruth4 needs to have image paths changed
%gTruth = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData_NoDebrisNOLIGHT_FULL_1944x2592\gTruth.mat');
gTruth = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData_NoDebrisNOLIGHT_FULL_1200x1600\gTruth.mat');
%gTruth = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData_DebrisNOLIGHT_FULL_1200x1600\gTruth.mat');
%gTruth = load('C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\sphereTrainingData_NoDebrisLIGHT_FULL_1200x1600\gTruth.mat'); %LIGHT
% pathPre = 'C:\Users\StevenSummey\Documents\MATLAB\Microscope Images\Testing_NO LIGHT';
% pathPost = [pwd '\Microscope Images\Training'];
% alternativePaths = [string(pathPre), string(pathPost)];
% unresolvedPaths = changeFilePaths(gTruth4.gTruth, alternativePaths);
% alternativePaths = [string(pathPost), string(pathPre)];

imds_1 = gTruth.gTruth.DataSource.Source;
imds = imageDatastore(imds_1);
imds.ReadFcn = @customReadDatastoreImage_NoDebris;
%imds.ReadFcn = @customReadDatastoreImage;

pxDir = gTruth.gTruth.LabelData.PixelLabelData;

classNames = ["Sphere", "Background"];
labelIDs = [1, 3];
pxds = pixelLabelDatastore(pxDir,classNames, labelIDs);
pxds.ReadFcn = @customReadDatastoreLabel_NoDebris;
%pxds.ReadFcn = @customReadDatastoreLabel;

targetSize = [1200 1600 3]; %smallest image dimension [972 1296 3], [486 648 3] [1944 2592 3], [1200 1600 3]

%%
%[imdsTrain, imdsVal, pxdsTrain, pxdsVal] = partitionImageData(imds_LIGHT,pxds_LIGHT); % LIGHT 
%[imdsTrain, imdsVal, pxdsTrain, pxdsVal] = partitionImageData(imds_NOLIGHT,pxds_NOLIGHT); % NO LIGHT
%[imdsTrain, imdsVal, pxdsTrain, pxdsVal] = partitionImageData_NoDebris(imds,pxds); % ALL, No Debris
[imdsTrain, imdsVal, pxdsTrain, pxdsVal] = partitionImageData_NoDebris(imds,pxds);
%[imdsTrain, imdsVal, pxdsTrain, pxdsVal] = partitionImageData(imds,pxds);
% [imdsTrain, imdsVal, imdsTest, pxdsTrain, pxdsVal, pxdsTest]
% Define validation data.
dsVal = combine(imdsVal,pxdsVal);
% Specify the network image size. This is typically the same as the traing image sizes.
imageSize = targetSize;

% Specify the number of classes.
numClasses = numel(classNames); 

% Create DeepLab v3+.
%lgraph18 = deeplabv3plusLayers(imageSize, numClasses, "resnet18");
%lgraph50 = deeplabv3plusLayers(imageSize, numClasses, "resnet50");
%lgraphMob = deeplabv3plusLayers(imageSize, numClasses, "mobilenetv2");
lgraphXcep = deeplabv3plusLayers(imageSize, numClasses, "xception");
%lgraphInc = deeplabv3plusLayers(imageSize, numClasses, "inceptionresnetv2");
% Define training options. 
options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropPeriod',5,...
    'LearnRateDropFactor',0.5,...
    'Momentum',0.9, ...
    'InitialLearnRate',0.01, ...
    'L2Regularization',0.005, ...
    'ValidationData',dsVal,...
    'MaxEpochs',30, ...  
    'MiniBatchSize',4, ...
    'Shuffle','every-epoch', ...
    'VerboseFrequency',2,...
    'Plots','training-progress',...
    'ValidationPatience', 4,...
    'CheckpointPath','C:\Users\StevenSummey\OneDrive - Aruna Bio\MATLAB\Image Processing\DeepLearningCheckpoint',...
    'CheckpointFrequency',5,...
    'CheckpointFrequencyUnit','epoch',...
    'ExecutionEnvironment','parallel');
%    'CheckpointPath', tempdir, ...
dsTrain = combine(imdsTrain, pxdsTrain);
xTrans = [-10 10];
yTrans = [-10 10];
dsTrain = transform(dsTrain, @(data)augmentImageAndLabel(data,xTrans,yTrans));
%%
%from checkpoint
net2 = trainNetwork(dsTrain, layerGraph(net), options);
%%
doTraining = true;
if doTraining    
    [net, info] = trainNetwork(dsTrain,lgraph,options);
end

%% TRAINING METRICS
    pxdsResults = semanticseg(imdsVal, trainedNetwork_1);

    metrics = evaluateSemanticSegmentation(pxdsResults, pxdsVal);
    
    figure;
    cm = confusionchart(metrics.ConfusionMatrix.Variables, ...
    classNames, Normalization='row-normalized');
    cm.Title = 'Normalized Confusion Matrix (%)';

    imageIoU = metrics.ImageMetrics.MeanIoU;
    figure
    histogram(imageIoU)
    title('Image Mean IoU')
%% TESTING
imDir = fullfile([pwd '\Microscope Images\Testing']);
imDirFiles = dir(fullfile(imDir,'*.tif'));
for i = 1:numel(imDirFiles)
    fullFile{i,1} = [imDirFiles(i).folder '\' imDirFiles(i).name];
end

imdsTEST = imageDatastore(fullFile); 
imdsTEST.ReadFcn = @customReadDatastoreImage;

for i = 1:numel(fullFile)
    figure;
    I = readimage(imdsTEST, i);
    C = semanticseg(I, trainedNetwork_1);
    BW = C == 'Sphere';
    B = labeloverlay(I,C,'Transparency',0.4);
    imshowpair(I,B,'montage')
end

%% 