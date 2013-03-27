function trainPatchModel(parDir, imNameList, selDim)

if( nargin < 1 )
    clear; close all; clc;
    display('No Input Arguments Provided');
    splitter = 1;
end
splitter = 0;
locRange = 1:18;
if( splitter == 1 )
 lbpRange = 19:58;
 grayRange = 59:130;
 hogRange = 131:194;
else
 lbpRange = 19:28;
 grayRange = 29:46;
 hogRange = 47:62;
end
if( nargin < 1 )
    parentDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/data/';
    posDirs = {'postSynapseDataset/synapseAlignedNew'};
    parDir = [parentDir posDirs{1}];
    load('fileNameStore');
    posShuffler = randperm( numPos );
    imNameList =  fileNames{1}(posShuffler(1:99));
    selDim = [locRange]; %194
end


posTrainCache = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/POS/'; system( ['rm ' posTrainCache '/*']);
negTrainCache = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/NEG/'; system( ['rm ' negTrainCache '/*']);
gtDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/groundTruth/';


imgIndex = [];
posFeatures = [];
negFeatures = [];
posLabel = [];
negLabel = [];
NUMRAND = 1;
patCtr = 1;
h = waitbar(0,'Please wait...');
for fileIter = 1:numel( imNameList )
    fileIter
    waitbar(fileIter/numel(imNameList),h);
    imName = imNameList(fileIter).name;
    fullFName{fileIter} = [parDir '/' imName];
    I = imread(fullFName{fileIter} );
    [noRows noCols noSlices] = size(I);
    patchSize = 256;
    xCenter = ceil( noCols/2 );
    yCenter = ceil( noRows/2 );
    xRange = [patchSize+1:noCols-patchSize]';
    yRange = [patchSize+1:noRows-patchSize]';
    [extension fName ext] = fileparts( imName );
    load([gtDir fName '.mat']);
    %gtPick  = [xCenter yCenter];
    gtPick = [ceil( (gtBox(1,1)+gtBox(1,3))/2 ) ceil( (gtBox(1,2)+gtBox(1,4))/2 )];
    randPick = [xRange( ceil( rand(NUMRAND,1)*numel(xRange) ) ) yRange( ceil( rand(NUMRAND,1)*numel(yRange) ) )];
    fullPicker = [gtPick; randPick];
    %imshow(I);
    %% Given a set of image with an assumption that the cetral patch comprises the target of interest,
    %% this code first samples all positive patches and stored them in a folder, along with random sampling of negative patches
    for imgIter = 1:NUMRAND+1
        imgIndex = [imgIndex; fileIter];
        % y1= min(1, (fullPicker(imgIter,2)-patchSize+1) );
        % y2= max(size(I,1), fullPicker(imgIter,2)+patchSize);
        % x1= min(1, (fullPicker(imgIter,1)-patchSize+1) );
        % x2= max(size(I,2), fullPicker(imgIter,1)+patchSize);
        y1= ( max(patchSize,fullPicker(imgIter,2)) -patchSize+1);
        y2= ( min(size(I,1)-patchSize, fullPicker(imgIter,2) ) +patchSize);
        x1= ( max(patchSize,fullPicker(imgIter,1))-patchSize+1);
        x2= ( min(size(I,2)-patchSize, fullPicker(imgIter,1) ) +patchSize);
        Icrop = I(y1:y2, x1:x2 );
        
        
        %% Features are extracted and stored for all these files
        localFeat = sd_cooccurFeature(Icrop, [], false);
        [lbpFeat, grayFeat, hogFeat] = compFeatures( Icrop, splitter );
        if( imgIter == 1 )
            %hold on; rectangle('Position', [x1 y1 x2-x1 y2-y1], 'LineWidth', 5, 'EdgeColor', [0 1 0] ); hold off;
            imTag = ['POS' num2str(imgIter)];
            fullPatchName = [posTrainCache fName imTag '.png'];
            imwrite( Icrop, fullPatchName );
            save(['/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/POS/' fName imTag '.mat'], 'lbpFeat', 'grayFeat', 'hogFeat', 'localFeat');
            posPatch{fileIter} = Icrop;
            posFeatures = [posFeatures; localFeat(:)' lbpFeat(:)' grayFeat(:)' hogFeat(:)'];
            posLabel = [posLabel; 1];
        else
            %hold on; rectangle('Position', [x1 y1 x2-x1 y2-y1], 'LineWidth', 5, 'EdgeColor', [1 0 0] ); hold off;
            imTag = ['NEG' num2str(imgIter)];
            negPatch{ (fileIter-1)*NUMRAND+imgIter } = Icrop;
            fullPatchName = [negTrainCache fName imTag '.png'];
            imwrite( Icrop, fullPatchName );
            save(['/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/NEG/' fName imTag '.mat'], 'lbpFeat', 'grayFeat', 'hogFeat', 'localFeat');
            negFeatures = [negFeatures; localFeat(:)' lbpFeat(:)' grayFeat(:)' hogFeat(:)'];
            negLabel = [negLabel; -1];
        end
        fullPName{patCtr} = fullPatchName;
        patCtr = patCtr + 1;
    end
end
close(h);
%selDim = [lbpRange grayRange hogRange];
fullLabels = [posLabel; negLabel];
fullFeatures = [posFeatures(:, selDim); negFeatures(:, selDim)];
save('testPatches', 'imgIndex', 'fullFeatures', 'fullLabels', 'fullPName');

%% SVM is learnt on different features
normMean = mean(fullFeatures);
normStd  = std(fullFeatures);
fullFeatures = ( fullFeatures - repmat( normMean, [size(fullFeatures,1) 1] ) )./repmat( normStd, [size(fullFeatures,1) 1] );
learntModel = train( fullLabels, sparse(fullFeatures) ); %, '-b 1 -t 2'

save('secondSynClassifier', 'learntModel', 'normMean', 'normStd', 'selDim', 'splitter');
[predLabels, decisionVals, predProbs] = predict( fullLabels, sparse(fullFeatures), learntModel, '-b 1' );
confusionmat( fullLabels, predLabels )