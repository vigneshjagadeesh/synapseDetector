clear; close all; clc;
addpath( genpath( 'C:\Users\vignesh\Dropbox\scratchSpace\synapseProject\interfaces\' ) );
if( ~ismac )
    patchStore = 'C:\researchCode\dataBase\synapseDataset\patchStore\';
else
    patchStore = '/media/OS/researchCode/dataBase/synapseDataset/patchStore/';
end
patchStore = '/media/OS/researchCode/dataBase/synapseDataset/patchStore/';
featVec = [];
featLabels = [];
VIS = true;
vesPatches = dir([patchStore 'Ves*']);
nonPatches = dir([patchStore 'Non*']);
fil = fspecial('log', 49, 7);
for iter = 1:numel(vesPatches)
    load( [patchStore vesPatches(iter).name] );
    if( ( size( c1, 1 ) + size( c1, 2 ) ) < 100 )
        continue;
    else
        [border1 mask1 ] = genBorderMask(c1 );
        [border2 mask2 ] = genBorderMask(c2 );
        c1 = ( adapthisteq(c1(:,:,1)) ); c2 = ( adapthisteq(c2(:,:,1)) );
        bw1 = im2bw(255-c1); bw2 = im2bw(255-c2); bw1 = imerode( bw1, strel('disk', 3) ); bw2 = imerode( bw2, strel('disk', 5) );
        
        
        log1 = imfilter( double(c1), fil, 'same', 'conv'); log2 = imfilter( double(c2), fil, 'same', 'conv');
        logMax1 = log1 > imdilate(log1, [1 1 1; 1 0 1; 1 1 1]); logMax2 = log2 > imdilate(log2, [1 1 1; 1 0 1; 1 1 1]);
        logMax1 = logMax1 .* border1; logMax2 = logMax2 .* border2;
        logFeat1(1) = sum( logMax1(:) ); logFeat2(1) = sum( logMax2(:) );
        for scaleIter = 1:4
            logFeat1( scaleIter + 1 ) = sum( sum( logMax1 .* mask1(:,:,scaleIter) .* log1 ) );
            logFeat2( scaleIter + 1 ) = sum( sum( logMax2 .* mask2(:,:,scaleIter) .* log2 ) );
        end
        if( VIS )
            figure(1); subplot(121); imshow(c1); subplot(122); imshow(c1); hold on; contour(bw1 , 'r' ); hold off;
            figure(2); subplot(121); imshow(c2); subplot(122); imshow(c2); hold on; contour(bw2 , 'r' ); hold off;
            figure(3); imagesc( log1 ); [y x] = find( logMax1 ); hold on; plot( x, y, 'g*'); hold off;
            figure(4); imagesc( log2 ); [y x] = find( logMax2 ); hold on; plot( x, y, 'g*'); hold off;
            pause(.25);
        end
        featVec = [featVec; logFeat1; logFeat2];
        featLabels = [featLabels; 1; 1];
    end
end

for iter = 1:numel(nonPatches)
    load( [patchStore nonPatches(iter).name] );
    if( ( size( c1, 1 ) + size( c1, 2 ) ) < 100 )
        continue;
    else
        [border1 mask1 ] = genBorderMask(c1 );
        [border2 mask2 ] = genBorderMask(c2 );
        c1 = ( adapthisteq(c1(:,:,1)) ); c2 = ( adapthisteq(c2(:,:,1)) );
        bw1 = im2bw(255-c1); bw2 = im2bw(255-c2); bw1 = imerode( bw1, strel('disk', 3) ); bw2 = imerode( bw2, strel('disk', 5) );
        
        log1 = imfilter( double(c1), fil, 'same', 'conv'); log2 = imfilter( double(c2), fil, 'same', 'conv');
        logMax1 = log1 > imdilate(log1, [1 1 1; 1 0 1; 1 1 1]); logMax2 = log2 > imdilate(log2, [1 1 1; 1 0 1; 1 1 1]);
        logMax1 = logMax1 .* border1; logMax2 = logMax2 .* border2;
        logFeat1(1) = sum( logMax1(:) ); logFeat2(1) = sum( logMax2(:) );
        for scaleIter = 1:4
            logFeat1( scaleIter + 1 ) = sum( sum( logMax1 .* mask1(:,:,scaleIter) .* log1 ) );
            logFeat2( scaleIter + 1 ) = sum( sum( logMax2 .* mask2(:,:,scaleIter) .* log2 ) );
        end
        if(VIS)
            figure(1); subplot(121); imshow(c1); subplot(122); imshow(c1); hold on; contour(bw1 , 'r' ); hold off;
            figure(2); subplot(121); imshow(c2); subplot(122); imshow(c2); hold on; contour(bw2 , 'r' ); hold off;
            figure(3); imagesc( log1 ); [y x] = find( logMax1 ); hold on; plot( x, y, 'g*'); hold off;
            figure(4); imagesc( log2 ); [y x] = find( logMax2 ); hold on; plot( x, y, 'g*'); hold off;
            pause(.25);
        end
        featVec = [featVec; logFeat1; logFeat2];
        featLabels = [featLabels; 2; 2];
    end
end
featVec = ( featVec - repmat( mean(featVec), size(featVec,1), 1 ) ) ./ repmat( std(featVec), size(featVec,1), 1 );
qualAssess = featureQuality( featVec, featLabels );






% Alg2
clear; close all; clc;
[trainDir trainFiles testDir testFiles FILTERINFO] = initializeParams();
USE_MIL = true;
%% Compute Training Features
trainCtr = 1;
classi = 1;
instLab = [];
for classIter = 1:2
    for fileIter = 3:10
        trainFileName{trainCtr} = [ trainDir{classIter} trainFiles{classIter}( fileIter ).name ];
        I = adapthisteq( rgb2gray( imread( trainFileName{trainCtr} ) ) );
        [imgbagTest trainBags{trainCtr, 1} X] = cbi_bagimg(I, FILTERINFO);
        if(classIter == 1)
            try,           labImg = imread([ trainDir{classIter} labFiles{classIter}( fileIter ).name ]);
            catch me,      labImg = zeros(size(I));
            end
            for iterX = 1:256:1024
                for iterY = 1:256:1024
                    instLab = [instLab; (unique( labImg(iterY:iterY+255, iterX:iterX+255) ) == 0) + 1];
                end
            end
        end
        trainBagLabels( trainCtr,1 ) = classIter;
        trainCtr = trainCtr + 1;
    end
    display([ 'Finished Class' num2str(classIter) ' in training features']);
end
%% Compute Testing Features
testCtr = 1;
for classIter = 1:2
    for fileIter = 3:30
        testFileName{testCtr} = [testDir{classIter} testFiles{classIter}(fileIter).name];
        try,
            I = adapthisteq( rgb2gray( imread( testFileName{testCtr} ) ) );
        catch me,
            display('PNG ERROR!');
        end
        [imgbagTest testBags{testCtr, 1} X] = cbi_bagimg(I, FILTERINFO);
        testBagLabels( testCtr,1 ) = classIter;
        testCtr = testCtr + 1;
    end
    display([ 'Finished Class' num2str(classIter) 'in testing features']);
end

% Whiten Data if Necessary
trainBags = whitenBagsSyn(trainBags, FILTERINFO.ID==1);
testBags = whitenBagsSyn(testBags, FILTERINFO.ID==1);

%% mi-Graph based multiple instance learning
if( USE_MIL )
    param.c = 100; param.gamma = .05; param.thresh = 40;
    trainKernelMat = constructKernel(trainBags, trainBags, param);
    testKernelMat = constructKernel(testBags, trainBags, param);
    opt = ['-c ',num2str(param.c),' -t 4 -s 0'];
    trainKernel = [(1:size(trainKernelMat, 1))' trainKernelMat];
    testKernel = [(1:size(testKernelMat, 1))' testKernelMat];
    model = svmtrain(trainBagLabels, trainKernel, opt);
    predLabels = svmpredict(testBagLabels, testKernel, model );
    figure(1);
    ctr = 1;
    for iter = 1:numel(testFileName), imshow(imread([testFileName{iter}])); title( [num2str(predLabels(iter) ) '  ' num2str(testBagLabels(iter) )]  ); pause(1); end
    pollerAcc = sum( predLabels(:) ~= testBagLabels(:) ) ./ numel(testBagLabels);
else
    USE_DIM_RED = false;
    %% Simple Instance Classifier with SVM,
    trainMat = cell2mat( trainBags );
    if( sum ( unique(instLab(:)) ) == 2 )
        trainMatLabel = [ones(sum(trainBagLabels==1)*size(trainBags{1}, 1), 1); 2*ones(sum(trainBagLabels==2)*size(trainBags{end}, 1), 1)];
    else
        trainMatLabel = [instLab(:); 2*ones(sum(trainBagLabels==2)*size(trainBags{end}, 1), 1)];
    end
    testMat  = cell2mat( testBags );
    testMatLabel = repmat( testBagLabels', size(testBags{1}, 1), 1);
    svmModel = svmtrain( trainMatLabel(:), trainMat );
    save( 'vesicle2dsvm19June11.mat', 'svmModel');
    if( classi == 1),      predMatLabel = svmpredict( testMatLabel(:), testMat, svmModel );
    elseif( classi == 2),  predMatLabel = double( boostClassifier( trainMat, trainMatLabel(:), testMat  ) );
    elseif( classi == 3 ), predMatLabel = classify( testMat, trainMat, trainMatLabel(:)  );
    else                   nbc = NaiveBayes.fit(trainMat,trainMatLabel(:)) ; predMatLabel = nbc.predict( testMat );
    end
    bagLabeler = reshape( predMatLabel, size(testBags{1}, 1), numel(testBags) );
    
    for iter = 1:size(bagLabeler, 2)
        currLab = reshape( bagLabeler(:,iter), [4 4] );
        currI =  adapthisteq( rgb2gray( imread(testFileName{iter}) ) );
        simBox(currI, currLab, 256); pause;
    end
    
    predLabels = max( bagLabeler, [], 1);
    pollerAcc = sum( predLabels(:) ~= testBagLabels(:) ) ./ numel(testBagLabels);
    % for iter = 1:numel(testFileName), createHeatMap(imread([testFileName{iter}]), bagLabeler(:,iter), FILTERINFO); pause; end
    if( USE_DIM_RED )
        %% Dimensionality Reduction for Visualization
        [coeff score] = princomp(trainMat);
        figure(1);
        posRed = score( trainMatLabel == 1, :);
        negRed = score( trainMatLabel == 2, :);
        plot( posRed(:,1), posRed(:,2), 'r*' );
        hold on;
        plot( negRed(:,1), negRed(:,2),  'g*' );
        hold off;
        
        figure(2);
        posRed = trainMat( trainMatLabel == 1, :);
        negRed = trainMat( trainMatLabel == 2, :);
        plot( posRed(:,1), posRed(:,2), 'r*' );
        hold on;
        plot( negRed(:,1), negRed(:,2),  'g*' );
        hold off;
        
        figure(3);
        for iter = 1:min(10, size(posRed,2) )
            subplot(5,2,iter);
            hold on;
            plot(posRed(:,iter), 'r*');
            plot(negRed(:,iter), 'g*');
            title( num2str(iter) );
            hold off;
        end
    end
    cvErr = dataQuality( trainMat, trainMatLabel);
end

% Alg 3
function logFeat1 = vesicleSearch(c1, fil)

%% Vesicle Features Heuristic threshold for current image size set to 1.5 top 500 retentions
% S1: Perform Median Filtering - Kernel 5 * 5
% S2: Adaptive Histogram Equalization
% S3: Perform Laplacian Gaussian Filtering
% S4: Threshold Weak Responses (set to 1.5)
% S5: Detect Peaks in 2D Map
% S6: Remove Border Artifacts
% S7: Retain top 500 responses by ranking them
% S8: Compute full Image Feature and Dispersion
% S9: Compute Pyramid Features of Dispersion

[border1 mask1 ] = genBorderMask(c1 );
masker = zeros( size(c1,1), size(c1,2) );
masker(40:end-40, 40:end-40 ) = 1;

c1 = medfilt2( double(c1), [5 5] );                                         %S1
cEq = double( adapthisteq( uint8(c1) ) );                                   %S2
log1 = imfilter( cEq, fil, 'same', 'conv');                                 %S3
log1( log1 < 1.5) = 0;                                                      %S4
logMax1 = log1 > imdilate(log1, [1 1 1; 1 0 1; 1 1 1]);                     %S5
logMax1 = logMax1 .* border1;                                               %S6
logMax1 = retainTop( log1, logMax1, 500);                                   %S7
[yF xF] = find( logMax1~=0 );                                               %S8
logFeat1(1) = sum( logMax1(:) .* log1(:) );
logDisp1(1) = median( pdist([xF yF],'euclidean') );
for scaleIter = 1:4                                                         %S9
    currMask = logMax1 .* mask1(:,:,scaleIter);
    [y x] = find( currMask~= 0 );
    logFeat1( scaleIter + 1 ) = sum( sum( currMask .* log1 ) );
    logDisp1( scaleIter + 1 ) = median( pdist([x y],'euclidean') );
end
logFeat1 = [logFeat1(:); logDisp1(:)];

function binMaskTop = retainTop( response, binMask, NUM)
binMaskTop = zeros( size(binMask,1), size(binMask,2) );
[y x]= find(binMask ~= 0 );
allVals = response( binMask ~= 0 );
[sortedRes sortedIds] = sort( allVals, 'descend' );
NUM = min( NUM, numel(y) );
for iter = 1:NUM
    binMaskTop( y( sortedIds(iter) ), x( sortedIds(iter) ) ) = 1;
end
