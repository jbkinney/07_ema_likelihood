clear all
close all

% complie all files
mex runmodelondata.cpp mi.cpp;
mex mifromscoresbybaseadjustment.cpp mi.cpp;
k = 1;

% Load data
dataLee = loadLeeData('abf1');
dataPerBin = 50;

% Remove 3/4 of data
N = numel(dataLee.regionIndices);
n = floor(N/16);
indices = randsample(N,n);
dataLee.regionScores = dataLee.regionScores(indices);
dataLee.regionIndices = dataLee.regionIndices(indices);


% abf1_lee_oneFourth
if true
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, numRuns);
    dataSets = repmat(quantizeData(dataLee, dataPerBin), 1, numRuns);
    computations(k) = ...
        MCMCComputation('abf1_lee_oneFourth', models, dataSets, dataLee); k=k+1;
end

% Choose runlocally based on whether or not FAFNER_HOME is set
if getenv('FAFNER_HOME')
    runlocally = false;
else
    runlocally = true;
end

% Que computations
quecomputations(computations, runlocally)