function [bw finalMask bbox localFeat genericFeat currImg] = simpleDetect(I, model, trTest, VIS, resStore);
%% Simple Function to Perform Synapse Localization
% Input Arguments:
% I - Input Image must be in multiples of 256
% model - Learnt Connected Component Model
% trTest - Either Perform Training orTesting (trTest = 1 means training)
% VIS - Give 1 for Visualization
% Use %alpha(h, finalMask256*.9); for visualizing an alpha map
% Output Arguments:
% bw - Binary Mask
% finalMask256 - Gives outputs in tile format
% bbox - Gives N detected bounding boxes N * [x1 y1 x2 y2]

% Info
%

%% If the user has provided no input arguments, try using an example
if( nargin < 1 )
    clear; close all; clc;
    model.minSize = 2667;
    model.maxSize = 15862;
    model.meanInt = 107.3966;
    model.stdInt  = 13.7951;
    trTest = 0;
    VIS = true;
    I = imread('/cluster/home/retinamap/datasetCreator/utahDataBase/postSynapseDataset/8_4_L274X181_Y161.png');
end
addpath('patchDescriptor');
[F FFreq] = fl_makeSubFilters(49); F = F(:,:,31:36);
fil = fspecial('log', 49, 7);
load('secondSynClassifier');
if( ~exist('splitter') ), splitter = 0; end

%% Intialize variables and the program, including output parameters
bbox = [];
genericFeat = [];
localFeat  = [];
currImg = [];
detScores = [];
[noRows noCols noSlices] = size(I);
finalMask = zeros(noRows, noCols );
R = 0;
thChoice = 2;
Iorig = I;
I = double(I);

%% Treshold Image to Create a Binary Mask
if( thChoice == 1 )   
    Ieq = histeq(Iorig);
    bw = Ieq < 20;
else
    Ieq = uint8( ceil( 255 .* ( I - min(I(:)) ) / max(I(:)) ) );
    mask = (I > 20 & I < 60 );
    bw = Ieq < mean(Ieq(mask==1)); 
    bw = Ieq < 50;
end

%% Perform Connected Components Analysis and remove trivial solutions and then dilate
CC = bwconncomp(bw);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = find(numPixels < 750); % This only removes trivial solutions
for iter = 1:numel(idx)
    bw(CC.PixelIdxList{idx(iter)}) = 0;
end
bw = imdilate( bw, strel('disk', 9 ) );
CClean = bwconncomp(bw);
numPixels = cellfun(@numel,CClean.PixelIdxList);

%% If testing mode, use learnt model to further prune connected compoenets not in a size range
if( trTest ==0 )
    % Have 4 * maxSize
    idx = find( (numPixels < model.minSize) | (numPixels > 4*model.maxSize) ); % This only removes trivial solutions
    for iter = 1:numel(idx)
        bw(CClean.PixelIdxList{idx(iter)}) = 0;
    end
    CClean = bwconncomp(bw);
elseif( VIS == 1 )
    f1 = figure(1); imshow(Ieq); hold on; contour( bw, [0 0], 'g' );  hold off;
    set( f1, 'Position', [2266           1        1135        1521] );
    return;
end

%% If testing mode,
if( CClean.NumObjects > 0 && CClean.NumObjects < 60 ) % was 35
    
    if(VIS)
        figure(2); h=imshow(uint8(I)); hold on; contour(bw, [0 0], 'r' ); hold off;
        title(['Number of Candidates are ' num2str(CClean.NumObjects)]);
    end
    %% Score the Bounding Boxes
    
    for objIter = 1:(CClean.NumObjects)
        currPixels = CClean.PixelIdxList{objIter};
        currScore(objIter) = mean( I(currPixels) );%normpdf( mean( I(currPixels) ), model.meanInt, model.stdInt );
        %currScore(objIter) = -log( currScore(objIter) );
        if( currScore(objIter) < 30 )
            currScore(objIter) = 200;
        end
    end
    
    %% Sort and Render Bounding Boxes based on score
    [currScores idx] = sort(currScore);
    featBox = [];
    
    objCtr = 0;
    for objIter = 1:(CClean.NumObjects)
        [y x] = ind2sub( [noRows noCols], CClean.PixelIdxList{idx(objIter)});
        xCen = ceil( mean(x) ); 
        yCen = ceil( mean(y) );
        patSize = 256;
        if( (xCen > patSize) & (xCen < noCols-patSize) & (yCen > patSize) & (yCen < noRows-patSize) )
            objCtr = objCtr + 1;
            currImg(:,:,1,objCtr) = I( (yCen-patSize+1):(yCen+patSize), (xCen-patSize+1):(xCen+patSize) );
            featBox = [featBox; (xCen-patSize+1) (yCen-patSize+1) (xCen+patSize) (yCen+patSize)];
            hold on; plot( xCen, yCen, 'r*'); hold off;
            croppedI = I((yCen-patSize+1):(yCen+patSize), (xCen-patSize+1):(xCen+patSize) );
            [lbpFeat, grayFeat, hogFeat] = compFeatures( croppedI, splitter );
            genericFeat = [genericFeat; lbpFeat(:)' grayFeat(:)' hogFeat(:)'];
        else
            continue;
        end
        if( objCtr < 16 )
            
            finalMask(( min(y)):(max(y)), (min(x)):(max(x)) ) = 1;
            bbox = [bbox; min(x) min(y) max(x) max(y)];
            detScores = [detScores; currScore(objIter)];
        else
            bbox = [bbox; min(x) min(y) max(x) max(y)];            
        end
    end
    if( objCtr > 0 ),
        if( selDim(1) > 1 )
            localFeat = sd_cooccurFeature(I, [], false, featBox);
        else
            localFeat = zeros( size(featBox,1), 18 ); 
        end
        fullFeatures = [localFeat genericFeat];
        fullFeatures = fullFeatures(:, selDim );
        fullFeatures = ( fullFeatures - repmat(normMean, [size(fullFeatures,1) 1]) )./repmat(normStd, [size(fullFeatures,1) 1]);
        [predLabels, decisionVals, predProbs] = predict( ones(size(fullFeatures,1), 1), sparse(fullFeatures), learntModel, '-b 1' );
        for boxIter = 1:numel(predLabels)
            if(VIS && (predLabels(boxIter) == 1) )
                drawRect( featBox(boxIter,1)+128, featBox(boxIter,2)+128, featBox(boxIter,3)-featBox(boxIter,1)-256, featBox(boxIter,4)-featBox(boxIter,2)-256, ['Detection ' num2str(predProbs(boxIter))], [0 1 0] );
            elseif( VIS && (predLabels(boxIter) ~= 1) )
                drawRect( featBox(boxIter,1)+128, featBox(boxIter,2)+128, featBox(boxIter,3)-featBox(boxIter,1)-256, featBox(boxIter,4)-featBox(boxIter,2)-256, ['Detection ' num2str(predProbs(boxIter))], [1 0 0] );
            end
        end
        save(resStore, 'bbox', 'featBox', 'predProbs', 'predLabels');
        %if(VIS) 
        %    montage( uint8(currImg) );
        %end
        bbox = bbox( predLabels == 1, : );
        featBox = featBox( predLabels == 1, :);
    end
end
