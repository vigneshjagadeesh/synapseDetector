clear; close all; clc;
addpath(genpath('/home/vignesh/Dropbox/spider/')); use_spider;
addpath('/home/vignesh/Dropbox/scratchSpace/synapseProject/interfaces/classifiersOffShelf/liblinear-1.8/');
% Test the feature subset selection with spider rouines ...
if( 0 )
    IMGSIZE = 256;
    bmark.trainDir{1}          = '/media/OS/researchCode/dataBase/synapseDataset/synapseDataset/ribbon/ribTrain/';%[bmark.dataDir  'synapseDetection/train/posSamples/']; %[bmark.dataDir 'labImg/'];
    bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'rib*']); bmark.trainFiles{1} = bmark.trainFiles{1}(3:50);
    bmark.trainDir{2}          = '/media/OS/researchCode/dataBase/synapseDataset/synapseDataset/random/randTrain/'; %[bmark.dataDir  'synapseDetection/train/negSamples/'];
    bmark.trainFiles{2}        = dir(bmark.trainDir{2});          bmark.trainFiles{2} = bmark.trainFiles{2}(3:50);
    bmark.testDir{1}           = '/media/OS/researchCode/dataBase/synapseDataset/synapseDataset/ribbon/ribTest/';%[bmark.dataDir 'synapseDetection/test/posSamples/'];
    bmark.testFiles{1}         = dir(bmark.testDir{1});  bmark.testFiles{1} = bmark.testFiles{1}(3:50);
    bmark.testDir{2}           = '/media/OS/researchCode/dataBase/synapseDataset/synapseDataset/random/randTest/'; %[bmark.dataDir 'synapseDetection/test/negSamples/'];
    bmark.testFiles{2}         = dir(bmark.testDir{2});  bmark.testFiles{2} = bmark.testFiles{2}(3:50);
    ctr =1;
    for filtSup = [27 49 63]
        F{ctr} = makeSubFilters(filtSup);
        F{ctr} = cat(3, F{ctr}(:,:,12:18), F{ctr}(:,:,31:36), F{ctr}(:,:,41:48));
        ctr = ctr + 1;
    end
    
    imgbag = [];
    featbag = [];
    trainLabel = [];
    %% Extract features on train files
    fullFeat = [];
    fullLabels = [];
    ctr = 1;
    for classIter = 1:2
        featurePool = [];
        for fileIter = 4:numel(bmark.trainFiles{classIter})
            display(['Processing File Number ' num2str(fileIter) ' in ' num2str(classIter) ' th class']);
            I = imread([bmark.trainDir{classIter} bmark.trainFiles{classIter}(fileIter).name]); I = double(imresize(I(:,:,1),[IMGSIZE IMGSIZE]));
            if( size(I,3) == 1 )
                I = cat(3, I, I, I);
            end
            [featVal featName] = featureBank(I,F);
            fullFeat = [fullFeat; featVal(:)'];
            fullLabels = [fullLabels; classIter];
        end
    end
    save('compFeatures', 'fullFeat', 'fullLabels');
end
load('compFeatures');
fullLabels( fullLabels == 2 ) = -1;
meanVec = mean(fullFeat); maxVec = std(fullFeat);
fullFeat = ( fullFeat - repmat(meanVec, [size(fullFeat,1) 1]) ) ./ repmat(maxVec, [size(fullFeat, 1) 1]);
%model = train(fullLabels, sparse(fullFeat), '-s 2');
%[predLabels accuracy probEst] = predict(fullLabels, sparse(fullFeat), model, '-b 1');
%cMat = confusionmat( predLabels, fullLabels);


d = data(fullFeat,fullLabels);
a = svm;
[tr a]=train(a,d);
display(['LOSS BEFORE SELECTION']); tr1=loss(tr,'class_loss')
CMat = confusionmat(tr.X, tr.Y);
accFull = sum(diag(CMat))/sum(CMat(:));
figure(1); plot(1:137, accFull, 'r', 'LineWidth',2); axis([0 200 0 2]);
ctr = 1;
for datIter = 10:10:size(fullFeat,2)
    a=rfe;
    a.feat=datIter;
    a.output_rank=0;
    [r,a]=train(a,d);
    display(['LOSS AFTER SELECTION']); r1=loss(r,'class_loss')
    CMat = confusionmat(r.X, r.Y);
    acc(ctr) = sum(diag(CMat))/sum(CMat(:));
    ctr = ctr + 1;
end
hold on; plot(10:10:size(fullFeat,2), acc, 'g', 'LineWidth',2);
sel = find(a.rank<100);
