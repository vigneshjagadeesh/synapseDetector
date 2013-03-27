clear; close all; clc;

%parentDir = '/cluster/home/retinamap/datasetCreator/utahDataBase/';
parentDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/data/';
posDirs = {'postSynapseDataset/synapseAlignedNew'};
negDirs = {'YAC_deepBorder', 'GAC_deepBorder', 'CBb4_deepBorder', 'CBb5_deepBorder', 'Gly+_deepBorder', 'MG_deepBorder',...
            'CBa2_deepBorder', 'CBab_deepBorder', 'HC_deepBorder', 'BC_deepBorder'};
resultsDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/detectionResults/';
numPosDirs = numel( posDirs );
numNegDirs = numel( negDirs );
filNames{1} = [];
filNames{2} = [];
numPos = 100;
numNeg = 100;
load('learntModel');
addpath('patchDescriptor');
addpath('/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/patchDescriptor');

system( ['rm ' resultsDir '*.mat'] );

if( exist( 'fileNameStore.mat', 'file' ) )
    load('fileNameStore');
else
    numPos = 0;
for posDirIter = 1:numel( posDirs )
    temp = dir( [parentDir posDirs{posDirIter} '/*.png'] );
    fileNames{1}( numPos+(1:numel(temp)) ) = temp;
    dirSample{1}( numPos+(1:numel(temp)) ) = posDirIter * ones( numel(temp,1), 1 );
    numPos = numPos + numel(temp);
end

numNeg = 0;
for negDirIter = 1:numel( negDirs )
    temp = dir( [parentDir negDirs{negDirIter} '/*.png'] );
    fileNames{2}( numNeg+(1:numel(temp)) ) = temp;
    dirSample{2}( numNeg+(1:numel(temp)) ) = negDirIter * ones( numel(temp,1), 1 );
    numNeg = numNeg + numel(temp);
end
save('fileNameStore', 'fileNames', 'dirSample', 'numPos', 'numNeg');
end

trTest = 0;
VIS = true;
posShuffler = randperm( numPos );
trainPatchModel([parentDir posDirs{ 1 }], fileNames{1}(posShuffler(1:99)), 19:62);
posShuffler = randperm( numPos );
for posIter = 1:50 %numel( posShuffler )
    posIter
    tic,
    currPick = posShuffler( posIter );
    I = imread( [parentDir posDirs{ dirSample{1}(currPick) } '/' fileNames{1}(currPick).name] );
    %I = imread( [parentDir negDirs{ dirSample{2}(currPick) } '/' fileNames{2}(currPick).name] );
    [bw{posIter} finalMask{posIter} bbox{posIter} localFeat{posIter} genericFeat{posIter} currImg{posIter}] = simpleDetect(I, model, trTest, VIS, ...
        [resultsDir fileNames{1}(currPick).name '_det.mat']);
    if( trTest==1 )
        f3 = figure(3); %set( f3, 'Position', [66           1        1135        1521] );
        [trainPatch{posIter} coor] = imcrop(I);
        coor = ceil( coor );
        trainMask{posIter} = detMask( coor(2):coor(2)+coor(4), coor(1):coor(1)+coor(3) );
    end
    toc
end

if( trTest==1 )
    for patIter = 1:numel(trainPatch)
        figure(4); subplot(211); imshow( trainPatch{patIter} );
        subplot(212); imshow( trainMask{patIter} ); %pause;
        CC = bwconncomp(trainMask{patIter});
        numPixels = cellfun(@numel,CC.PixelIdxList);
        accStats(patIter,:) = [ mean( trainPatch{patIter}(:) ) max( numPixels ) ];
    end
    model.minSize = min( accStats(:,2 ) );
    model.maxSize = max( accStats(:,2 ) );
    model.meanInt = mean( accStats(:,1) );
    model.stdInt  = std( accStats(:,1) );
    save('learntModel', 'model');
    save('trainPatches', 'trainPatch', 'trainMask');
end

apDetections(resultsDir);