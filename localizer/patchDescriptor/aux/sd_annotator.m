clear; close all; clc;
bmark.parDir = 'C:/Users/vignesh/Dropbox/synapseProject/';        
bmark.dataDir = 'C:/researchCode/dataBase/synapseDataset/';
currDir = pwd; cd(bmark.dataDir); mkdir('patchStore'); cd( currDir );
patchStore = [bmark.dataDir 'patchStore/'];


% % % % % % bmark.trainDir{1}          = [bmark.dataDir  'synapseDataset/conventional/conSynS1/']; %[bmark.dataDir 'labImg/'];
% % % % % % bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'con*']); shuffler = randperm(numel(bmark.trainFiles{1}));
% % % % % % bmark.trainFiles{1} = bmark.trainFiles{1}( shuffler(1:50) );
% % % % % % figure(1);
% % % % % % for iter = 1:numel(bmark.trainFiles{1})
% % % % % %     currFileName = [bmark.trainDir{1} bmark.trainFiles{1}(iter).name];
% % % % % %     I = imread( currFileName ); 
% % % % % %     [c1 rect1] = imcrop(I); title('Annotate 1 vesicular patches');
% % % % % %     [c2 rect2] = imcrop(I); title('Annotate 2 vesicular patches');
% % % % % %     save([patchStore 'VesicularPatch ' num2str(iter) '.mat'], 'currFileName', 'c1', 'c2', 'rect1', 'rect2');
% % % % % % end

bmark.trainDir{1}          = [bmark.dataDir  'synapseDataset/random/randPatS1/']; %[bmark.dataDir 'labImg/'];
bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'rand*']); shuffler = randperm(numel(bmark.trainFiles{1}));
bmark.trainFiles{1} = bmark.trainFiles{1}(shuffler(1:50));
for iter = 1:numel(bmark.trainFiles{1})
    currFileName = [bmark.trainDir{1} bmark.trainFiles{1}(iter).name];
    I = imread( currFileName ); 
    [c1 rect1]  = imcrop(I); title('Annotate 1 non-vesicular patches');
    [c2 rect2] = imcrop(I); title('Annotate 2 non-vesicular patches');
    save([patchStore 'NonVesicularPatch ' num2str(iter) '.mat'], 'currFileName', 'c1', 'c2', 'rect1', 'rect2');
end
pause;
bmark.trainDir{1}          = [bmark.dataDir  'synapseDataset/ribbon/ribSynS1/']; %[bmark.dataDir 'labImg/'];
bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'con*']); shuffler = randperm(numel(bmark.trainFiles{1}));
bmark.trainFiles{1} = bmark.trainFiles{1}(shuffler(1:50));
for iter = 1:numel(bmark.trainFiles{1})
    currFileName = bmark.trainFiles{1}(iter).name;
    I = imread( currFileName ); I = histeq( I(:,:,1) );
    [c1 rect1]  = imcrop(I); title('Annotate 1 ribbon patches');
    [c2 rect2] = imcrop(I); title('Annotate 2 ribbon patches');
    save([patchStore 'RibbonPatch ' num2str(iter) '.mat'], 'currFileName', 'c1', 'c2');
end

bmark.trainDir{1}          = [bmark.dataDir  'synapseDataset/conventional/conSynS1/']; %[bmark.dataDir 'labImg/'];
bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'con*']); shuffler = randperm(numel(bmark.trainFiles{1}));
bmark.trainFiles{1} = bmark.trainFiles{1}(shuffler(1:50));
for iter = 1:numel(bmark.trainFiles{1})
    currFileName = bmark.trainFiles{1}(iter).name;
    I = imread( currFileName ); I = histeq( I(:,:,1) );
    [c1 rect1]  = imcrop(I); title('Annotate 1 cleft patches');
    [c2 rect2] = imcrop(I); title('Annotate 2 cleft patches');
    save([patchStore 'CleftPatch ' num2str(iter) '.mat'], 'currFileName', 'c1', 'c2');
end

bmark.trainDir{1}          = [bmark.dataDir  'synapseDataset/random/randPatS1/']; %[bmark.dataDir 'labImg/'];
bmark.trainFiles{1}        = dir([bmark.trainDir{1} 'con*']); shuffler = randperm(numel(bmark.trainFiles{1}));
bmark.trainFiles{1} = bmark.trainFiles{1}(shuffler(1:50));
for iter = 1:numel(bmark.trainFiles{1})
    currFileName = bmark.trainFiles{1}(iter).name;
    I = imread( currFileName ); I = histeq( I(:,:,1) );
    [c1 rect1]  = imcrop(I); title('Annotate 1 non-cleft patches');
    [c2 rect2] = imcrop(I); title('Annotate 2 non-cleft patches');
    save([patchStore 'NonCleftPatch ' num2str(iter) '.mat'], 'currFileName', 'c1', 'c2');
end

