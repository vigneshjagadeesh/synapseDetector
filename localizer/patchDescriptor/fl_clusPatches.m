%function predLabels = clusPatches(initNumClus, treeDepth, fullImg)
clear; close all; clc;
load('testPatches');
%save('testPatches', 'imgIndex', 'fullFeatures', 'fullLabels', 'fullFName');

%if( nargin < 1 )
    initNumClus = 6;
    treeDepth = 4;
%end

CHROME_VIS = true;
if( CHROME_VIS )
    system(['pkill chrome']);
    system(['/usr/bin/google-chrome http://www.google.com --new-window &']);
    display('Opened a Window of Chrome');
end

% Perform unsupervised clustering first
currDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/';
% % fullFeat = [];
% % fullName = [];
% % fullCoor = [];
% % fullLabel = [];
% % numImg = numel(featbag);
% % ctr = 1;
% % for imgIter = 1:numImg
% %     fullFeat = [fullFeat; featbag{imgIter}.feats];
% %     numPatches = size(featbag{imgIter}.feats,1);
% %     featbag{imgIter}.XFor = [featbag{imgIter}.XFor imgIter*ones(numPatches,1)];
% %     fullCoor = [fullCoor; featbag{imgIter}.XFor];
% %     for patIter = 1:numPatches
% %         fullName{ctr} = featbag{imgIter}.imgName;
% %         fullLabel = [fullLabel; featbag{imgIter}.trainLabel];
% %         ctr = ctr + 1;
% %     end
% % end
fullCoor = imgIndex;
fullName = fullPName;
trainLabel = fullLabels;
fullFeat = fullFeatures;
meanVec = mean(fullFeat);
maxVec = max(fullFeat);
%fullFeat = ( fullFeat - repmat(meanVec, [size(fullFeat,1) 1]) ) ./ repmat(maxVec, [size(fullFeat, 1) 1]);
trainLabel(trainLabel == 2) = -1;
trainLabelOrig = trainLabel;
fullDataOrig = fullFeat;

% Perform KMeans Clustering
latentTracker = ones( numel(trainLabel), 1 );
survivors = 1:numel(latentTracker);
for depthIter = 1:treeDepth
    numClus = ceil(initNumClus/depthIter); %max(3, ceil(initNumClus / depthIter) );
    system('rm -r patchCluster*');
    display(['Processing Depth ' num2str(depthIter) ' of tree']);
    currFeat = fullFeat(latentTracker==1, :);
    currLabels = trainLabel(latentTracker==1, :);
    origMapper = find( latentTracker == 1 );
    if(    ( numel(intersect(origMapper, survivors)) ~= numel(survivors) ) | ( numel(survivors) ~= numel(origMapper) )    )
        error('Elements do not match up');
    else
        display(['Currently having ' num2str(numel(survivors)) ' survivors']);
    end
    display(['Number of datapoints used in clustering is ' num2str(numel(origMapper)) ]);
    [idx C] = kmeans( currFeat, numClus,  'Replicates', 10);
    display(['Finished Clustering in Depth ' num2str(depthIter) ' of tree']);
    posProp = zeros(numClus, 1);
    removeIndices = [];
    for clusIter = 1:numClus
        currSelect = idx == clusIter;
        clusLab = currLabels( currSelect );
        posProp(clusIter) = sum(clusLab == 1)/numel(clusLab);
    end
    [sortVal sortInd] = sort( posProp, 'descend' );
    removeCand = sortInd( ceil(numel(sortInd)/2):end );
    for remIter = 1:numel(removeCand)
        removeIndices = [removeIndices; find(idx==removeCand(remIter))];
    end
    latentTracker( origMapper(removeIndices) ) = -1; % Fix up the latent tracker
    %latentTracker( currLabels == -1 ) = -1;
    display(['Estimated Postive Label Propotion' num2str(sum(latentTracker==1)/numel(latentTracker))] );
    scoreProgress(depthIter) = mean(posProp);
    
    % Visualization Rotuine
    survivors = [];
    ctr = 1;
    for iter = 1:numClus
        if( sum(removeCand == iter) > 0 ), clusterStatus = 'Reject'; else, clusterStatus = 'Accept'; end
        mkdir(['patchCluster' num2str(iter)]);
        currInd = find( idx == iter ); % This is all patches in cluster 1
        for patIter = 1:numel(currInd)
            currIndex = origMapper( currInd(patIter) );
            if( strcmp(clusterStatus, 'Accept') )
                survivors = [survivors; currIndex];
            end
            display(['Patch ' num2str(ctr) ' has index ' num2str(currIndex) ' and is ' clusterStatus]); ctr = ctr + 1;            
            currPatch = imread(fullPName{ceil(currIndex)});
            %currPatch = fullImg{fullCoor(currIndex,5)}(fullCoor(currIndex,2):fullCoor(currIndex,2)+fullCoor(currIndex,4)-1, fullCoor(currIndex,1):...
            %                    fullCoor(currIndex,1)+fullCoor(currIndex,3)-1, 1);
            currPatch(1:64, 1:64, : ) = 255 * (trainLabelOrig(currIndex)==1);
            imwrite(uint8(currPatch), ['patchCluster' num2str(iter) '/patch' num2str(patIter) '.jpg'] );
        end
    end
    
    if( CHROME_VIS )
        for iter = 1:numClus
            if( sum(removeCand == iter) > 0 ), clusterStatus = 'Reject'; else, clusterStatus = 'Accept'; end
            system(['python createHtml.py ' currDir 'patchCluster' num2str(iter) '/ Depth:' num2str(depthIter) '--Cluster_Num:' num2str(iter) '--Cluster_Status:' clusterStatus ]);
            if( ismac )
                system(['/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome ' currDir 'patchCluster' num2str(iter) 'sampleViewer.html']);
            else
                %system(['/usr/bin/google-chrome ' currDir 'patchCluster' num2str(iter) '/sampleViewer.html']);
                web([currDir 'patchCluster' num2str(iter) '/sampleViewer.html'], '-new');
            end
        end
        %pause;
    end
end
model = train(latentTracker, sparse(fullFeat), '-s 2');
[predLabels accuracy probEst] = predict(latentTracker, sparse(fullFeat), model, '-b 1');