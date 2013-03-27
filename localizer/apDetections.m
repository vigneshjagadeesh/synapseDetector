%function apDetections(resultsDir)
clear; close all; clc;
resultsDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/detectionResults/';
%
if( nargin < 1 )
    resultsDir = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/detectionResults/';
end
allFiles = dir([resultsDir '*.mat']);
posDirs = {'/cluster/home/vignesh/retinamapClone/multiAttributes/data/postSynapseDataset/synapseAlignedNew/'};
gtDirs  = '/cluster/home/vignesh/retinamapClone/multiAttributes/localizer/trainCache/groundTruth/';
algDetections = [];
gtDetections = [];
algScores = [];
algIndicators = [];
tp = [];
fp = [];
numpos = 0;
for fileIter = 1:numel( allFiles )
    flagger = 1;
    [fileName fileExt] = strtok( allFiles(fileIter).name, 'png' );
    canvasImageGT = zeros(2559, 2559 );
    gtDetections = [gtDetections; 1023 1023 1535 1535];
    canvasImageGT(1023:1535, 1023:1535 ) = 1;
    %gtDetections = [gtDetections; 1151 1151 1407 1407];
    %canvasImageGT(1151:1407, 1151:1407 ) = 1;
    load( [resultsDir allFiles(fileIter).name] );
    load( [gtDirs fileName 'mat'] );
    I = imread( [posDirs{1} fileName 'png'] );
    imshow(I);
    for iter = 1:size(bbox,1)
        hold on; drawRect(bbox(iter,1), bbox(iter,2), bbox(iter,3)-bbox(iter,1), bbox(iter,4)-bbox(iter,2), '1', [1 0 0] ); hold off;
    end
    for iter = 1:size(gtBox,1)
        hold on; drawRect(gtBox(iter,1), gtBox(iter,2), gtBox(iter,3)-gtBox(iter,1), gtBox(iter,4)-gtBox(iter,2), '1', [0 1 0] ); hold off;
    end
    
    try,
        featBox = featBox( predLabels == 1, : );
        predProbs = predProbs(predLabels == 1,1);
    catch me,
        continue;
    end
    numpos = numpos + 1;
    algDetections = [algDetections; featBox];
    algScores = [algScores; predProbs];
    bbgt = gtBox(1, :);
    for boxIter = 1:size(featBox,1)
        bb = featBox(boxIter, :);
        canvasImageAL = zeros(2559, 2559 );
        canvasImageAL(bb(2):bb(4), bb(1):bb(3)) = 1;
        tempMask = canvasImageAL(:) + canvasImageGT(:);
        overlapScore = sum( tempMask(:) == 2 )/sum(tempMask(:)>0);
        
        bi=[max(bb(1),bbgt(1)) ; max(bb(2),bbgt(2)) ; min(bb(3),bbgt(3)) ; min(bb(4),bbgt(4))];
        iw=bi(3)-bi(1)+1;
        ih=bi(4)-bi(2)+1;
        if iw>0 & ih>0
            % compute overlap as area of intersection / area of union
            ua=(bb(3)-bb(1)+1)*(bb(4)-bb(2)+1)+...
                (bbgt(3)-bbgt(1)+1)*(bbgt(4)-bbgt(2)+1)-...
                iw*ih;
            overlapScore=iw*ih/ua;
        else
            overlapScore = 0;
        end
        
        
        algIndicators = [algIndicators; (overlapScore>0)];
        if( (algIndicators(end) == 1) && (flagger == 1) )
            tp = [tp;1];
            fp =[fp; 0];
            flagger = 0;
        else
            if( rand > 0 )
                tp = [tp;0];
                fp = [fp;1];
            else
                numpos = numpos + 1;
                tp = [tp;1];
                fp =[fp; 0];
            end
        end
    end
    
end

[sortVal sortIdx] = sort( algScores );
tpSorted = tp( sortIdx);
fpSorted = fp( sortIdx);
fp=cumsum(fpSorted);
tp=cumsum(tpSorted);
rec=tp/numpos;
prec=tp./(fp+tp);
figure(1); plot( rec, prec, 'r' ); grid on;
mrec=[0 ; rec ; 1];
mpre=[0 ; prec ; 0];
for i=numel(mpre)-1:-1:1
    mpre(i)=max(mpre(i),mpre(i+1));
end
i=find(mrec(2:end)~=mrec(1:end-1))+1;
ap=sum((mrec(i)-mrec(i-1)).*mpre(i))
