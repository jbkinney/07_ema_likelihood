clear all
close all

% complie all files
mex runmodelondata.cpp mi.cpp;
mex mifromscoresbybaseadjustment.cpp mi.cpp;
k = 1;

% Load data
dataMuk = loadMukData('abf1');
dataLee = loadLeeData('abf1');
dataPerBin = 50;

% abf1_muk
if true
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, numRuns);
    dataSets = repmat(quantizeData(dataMuk, dataPerBin), 1, numRuns);
    computations(k) = ...
        MCMCComputation('abf1_muk', models, dataSets, dataMuk); k=k+1;
end

% abf1_lee
if false
    clear('dataSets', 'models');
    numRuns = 10;
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, numRuns);
    dataSets = repmat(quantizeData(dataLee, dataPerBin), 1, numRuns);
    computations(k) = ...
        MCMCComputation('abf1_lee', models, dataSets, dataLee); k=k+1;
end

% abf1_lee_halfA and abf1_lee_halfB
if false
    clear('dataSetsA', 'dataSetsB', 'models');
    numRuns = 5;
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, numRuns);
    [dataSetA, dataSetB] = partitionData(dataLee, dataPerBin);
    dataSetsA = repmat(dataSetA,1,5);
    dataSetsB = repmat(dataSetB,1,5);   
    computations(k) = ...
        MCMCComputation('abf1_lee_halfA', models, dataSetsA, dataLee); k=k+1;
    computations(k) = ...
        MCMCComputation('abf1_lee_halfB', models, dataSetsB, dataLee); k=k+1;
end

% abf1_lee_widths
if false
    clear('dataSets', 'models');
    numRuns = 1;
    motifs = {
       'NNNNNNRTCAYTNNNNACGWNNNNNN';
        'NNNNNRTCAYTNNNNACGWNNNNN'; 
         'NNNNRTCAYTNNNNACGWNNNN';
          'NNNRTCAYTNNNNACGWNNN';
           'NNRTCAYTNNNNACGWNN';
            'NRTCAYTNNNNACGWN';
             'RTCAYTNNNNACGW';
              'TCAYTNNNNACG';
               'CAYTNNNNAC';
                'AYTNNNNA'};  
    N = numel(motifs);
    for i=1:N
        models(i) = motif2model(motifs{i});
    end
    dataSets = repmat(quantizeData(dataLee, dataPerBin), 1, N);
    computations(k) = ...
        MCMCComputation('abf1_lee_widths', models, dataSets, dataLee); k=k+1;
end

% abf1_lee_bins
if false
    clear('dataSets', 'models');
    numRuns = 1;
    dataPerBinChoices = [20, 30, 40, 50, 75, 100, 125, 150, 200];
    N = numel(dataPerBinChoices);
    for i=1:N
        dataSets(1,i) = quantizeData(dataLee, dataPerBinChoices(i));
    end
    models = repmat(motif2model('NNNRTCAYTNNNNACGWNNN'), 1, N);
    computations(k) = ...
        MCMCComputation('abf1_lee_bins', models, dataSets, dataLee); k=k+1;
end

% Choose runlocally based on whether or not FAFNER_HOME is set
if getenv('FAFNER_HOME')
    runlocally = false;
else
    runlocally = true;
end

% Que computations
quecomputations(computations, runlocally)