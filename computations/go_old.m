clear all
close all

% complie all files
mex runmodelondata.cpp mi.cpp;
mex mifromscoresbybaseadjustment.cpp mi.cpp;
k = 1;

% rap1_muk_20
if false
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNWRMACCCATACAYYNNN'), 1, numRuns);
    data = loadMukData('rap1');
    dataSets = repmat(quantizeData(data, 20), 1, numRuns);
    computations(k) = ...
        MCMCComputation('rap1_muk_20', models, dataSets, data); k=k+1;
end

% rap1_lee_20
if false
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNWRMACCCATACAYYNNN'), 1, numRuns);
    data = loadLeeData('rap1');
    dataSets = repmat(quantizeData(data, 20), 1, numRuns);
    computations(k) = ...
        MCMCComputation('rap1_lee_20', models, dataSets, data); k=k+1;
end

% reb1_lee_50
if true
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNNTTACCCGNNNN'), 1, numRuns);
    data = loadLeeData('reb1');
    dataSets = repmat(quantizeData(data, 50), 1, numRuns);
    computations(k) = ...
        MCMCComputation('reb1_lee_50', models, dataSets, data); k=k+1;
end

% reb1_lee_bins
if true
    clear('dataSets', 'models');
    numRuns = 1;
    dataPerBinChoices = [20, 30, 40, 50, 75, 100, 125, 150, 200];
    N = numel(dataPerBinChoices);
    data = loadLeeData('reb1');
    for i=1:N
        dataSets(1,i) = quantizeData(data, dataPerBinChoices(i));
    end
    models = repmat(motif2model('NNNNTTACCCGNNNN'), 1, N);
    computations(k) = ...
        MCMCComputation('reb1_lee_bins', models, dataSets, data); k=k+1;
end

% cin5_lee_50
if false
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNNTTACRTAANNN'), 1, numRuns);
    data = loadLeeData('cin5');
    dataSets = repmat(quantizeData(data, 50), 1, numRuns);
    computations(k) = ...
        MCMCComputation('cin5_lee_50', models, dataSets, data); k=k+1;
end


% Choose runlocally based on whether or not FAFNER_HOME is set
if getenv('FAFNER_HOME')
    runlocally = false;
else
    runlocally = true;
end

% Que computations
quecomputations(computations, runlocally)