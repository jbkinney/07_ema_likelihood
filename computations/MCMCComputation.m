% A seed consists of
% model0, dataSet, data
% Multiple dataSets should be aligned in ROWS!
function c = MCMCComputation(name, models, dataSets, data)
    
    % Make sure each model has its own data set (aligned along rows)
    if size(dataSets,2) ~= numel(models)
        error('size(dataSets,2) ~= numel(models)')
    end

    c.description = name;
    c.name = name; 
    c.dir = ['mcmc_results/' c.name];
    c.laciDir = '~jbkinney/Research/EMALikelihood/computations/mcmc_results';
    c.recipients = 'jbkinney@princeton.edu';
    c.data = data;
    c.numDataSets = numel(c.data);
    c.numRuns = numel(models); 
    
    for m=1:c.numRuns
        r.fileName = ['run_' num2str(m) '.mat'];
        r.finishedFile = ['finished_' num2str(m) '.mat'];
        
        r.model0 = fixModelGauge(models(m));
        r.dataSet = dataSets(:,m);
        r.runNum = m;
        r.L = size(r.model0.emat,1);
        r.maxNumTrials = 2e6; %5e6;
        r.trialsPerWrite = 10000;
        r.trialsPerJot = 1000;
        r.trialsPerStep = 100;
        c.runs(m,1) = r;
    end
    
    % Set main computation funciton
    c.main = @(runNum, runlocally) mcmc(c,runNum,runlocally);
end