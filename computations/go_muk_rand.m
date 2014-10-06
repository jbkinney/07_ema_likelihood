clear all
close all

% complie all files
mex runmodelondata.cpp mi.cpp;
mex mifromscoresbybaseadjustment.cpp mi.cpp;
k = 1;

% Load data
dataMuk = loadMukData('abf1');
dataPerBin = 20;

% Randomize data
scores = dataMuk.regionScores;
dataMuk.regionScores = randsample(scores, numel(scores));

% abf1_muk_rand_20
if true
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, numRuns);
    dataSets = repmat(quantizeData(dataMuk, dataPerBin), 1, numRuns);
    computations(k) = ...
        MCMCComputation('abf1_muk_rand_20', models, dataSets, dataMuk); k=k+1;
end

% Choose runlocally based on whether or not FAFNER_HOME is set
if getenv('FAFNER_HOME')
    runlocally = false;
else
    runlocally = true;
end

% Que computations
quecomputations(computations, runlocally)