function [trainDir trainFiles testDir testFiles FILTERINFO] = initializeParams()
if( ispc )
    parDir = 'C:/Users/vignesh/Dropbox/synapseDetect/';
    dataDir = 'C:/researchCode/dataBase/synapseDataset/';
else
    parDir = '/Users/vigneshjagadeesh/Dropbox/synapseDetect/';
    dataDir = '/Users/vigneshjagadeesh/Dropbox/synapseDetect/synapseDataset/';
end
addpath( genpath( [parDir 'interfaces'] ) );
trainDir{1} = [dataDir 'controlSet/train/posSamples/'];       trainFiles{1} = dir([trainDir{1} 'con*']); labFiles{1} = dir([trainDir{1} 'Label*']);
trainDir{2} = [dataDir  'controlSet/train/negSamples/'];      trainFiles{2} = dir(trainDir{2});
testDir{1}  =  [dataDir 'controlSet/test/posSamples/']; testFiles{1} = dir(testDir{1});
testDir{2}  =  [dataDir 'controlSet/test/negSamples/']; testFiles{2} = dir(testDir{2});
FILTERINFO.N = 256;
FILTERINFO.F = cbi_gabordictionary(6, 6, FILTERINFO.N, [0.05 0.4], 0); %
FILTERINFO.gistFilt = createGabor([8 8 4], FILTERINFO.N);
FILTERINFO.ID = 1;
FILTERINFO.numberGistBlocks = 4;
FILTERINFO.patchType = 'nonoverlap';
FILTERINFO.mapping=getmapping(8,'u2');