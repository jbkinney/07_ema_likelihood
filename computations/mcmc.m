% Main MCMC routine
function mcmc(c, runNum, runlocally)
    %rand('seed', runNum);   % MUST SEED RANDOM NUM GENERATOR!!!

    if nargin == 2
        runlocally = false;
    end

    if runlocally
        figure('windowstyle', 'docked');
    end
    
    r = c.runs(runNum);
    numJots = floor(r.maxNumTrials/r.trialsPerJot);
    for i=1:numel(r.dataSet)
        r.N(i,1) = numel(r.dataSet(i).regionTypes);
    end
    Is = zeros(numJots, c.numDataSets);
    model = fixModelGauge(r.model0);    
    models = repmat(model, [numJots, 1]);
    allIs = zeros(r.maxNumTrials, c.numDataSets);
  
    % Make sure initial model is within bounds
    oldI = runmodelondata(model, r.dataSet, c.data, 0);
    I = oldI;

    % Run MCMC loop for r.maxNumTrials
    trialNum = 0;
    jotNum = 0;
    r.numRejects = 0;
    while (trialNum < r.maxNumTrials)

        % Choose position if not specified by user
        x = randsample(r.L,2);
        pos1 = x(1);
        pos2 = x(2);

        % Set scoresByBase
        [I, regionStats] = runmodelondata(model, r.dataSet, c.data, 0, pos1, pos2);
        
        if trialNum < r.maxNumTrials/4
            sigma = .1;
        else
            sigma = .02;
        end
        
        
        % Turbo MCMC
        adjustments = zeros(16,1);
        numTrialsInThisStep = poissrnd(r.trialsPerStep);
        for t=1:numTrialsInThisStep 
            trialNum = trialNum+1;
            newmodel = model;
            
            % If adjusting an emat element...
            reject = false;
            randNum  = rand();
            if randNum < .8
            
                % To adjust an emat element, first chose an i and b 
                i = randsample([pos1, pos2],1); %Choose i = pos1 or pos2
                b = randsample(1:4,1);          % Choose base
            
                % Get correpsnding cols
                if i == pos1
                    cols = 4*(b-1) + (1:4);
                else
                    cols = 4*((1:4)-1) + b;
                end
                
                % Adjust the energy matrix. Reject if out of range
                d = sigma*randn();
                y = newmodel.emat(i,b)+d;
                newmodel.emat(i,b) = y;
                if (y > 1) | (y < 0)
                    reject = true;
                end

            % Adjusting a cutoff
            elseif randNum < .9
                % Choose a cutoff
                n = randsample(c.numDataSets,1);
                cut = newmodel.cutoff(n);
                newcut = cut + sigma*randn();
                newmodel.cutoff(n) = newcut;  
                
                cols = 1:16;
                d = 0;
            % Otherwise, adjust entire column and cutoff
            else
                cols = 1:16;
                d = sigma*randn();
                for n=1:c.numDataSets;
                    newmodel.cutoff(n) = newmodel.cutoff(n)+d;
                end
                i = randsample([pos1, pos2], 1);
                row = newmodel.emat(i,:)+d;
                newmodel.emat(i,:) = row;
                if (max(row) > 1) | (min(row) < 0)
                    reject = true;
                end
            end
        
            % Record adjustments for quick recomputation
            adjustments(cols) = adjustments(cols)+d;
            
            % Only compute I if not rejecting
            if ~reject
                %I = runmodelondata(newmodel, r.dataSet, c.data, 0, pos1, pos2);
                for n=1:numel(r.dataSet)  
                    I(n) = mifromscoresbybaseadjustment(adjustments, ...
                        regionStats(n).scoresByBase, newmodel.cutoff(n), ...
                        r.dataSet(n).regionTypes, r.dataSet(n).numTypes);
                end
            end
            
            % Decide whether or not to keep new mutual info
            p = exp(r.N.*(I-oldI));
            if (~reject) & (rand() < p)
                model = newmodel;
                oldI = I;
            else
                r.numRejects = r.numRejects + 1;
                I = oldI;
                adjustments(cols) = adjustments(cols)-d;
            end
            allIs(trialNum,:) = I;
            
            % Record model and Fs.
            if mod(trialNum,r.trialsPerJot) == 0
                jotNum = jotNum+1;
                models(jotNum) = model;
                Is(jotNum,:) = I(:);
                
                if runlocally
                    rejectFraction = r.numRejects/trialNum
                    subplot(4,1,1)
                    gaugemodel = fixModelGauge(model);
                    imagesc(gaugemodel.emat');
                    title(['cutoff == ' num2str(gaugemodel.cutoff') ', trialNum == ' num2str(trialNum)]);
                    colorbar;
            
                    subplot(4,1,2)
                    z = allIs(1:trialNum, 1);
                    %z = z(z >= (max(z) - .02));
                    %[y, x] = hist(z, 50);
                    minz = min(max(allIs)'-.02);
                    %maxz = max(allIs(:));
                    [y,x] = hist(allIs(min(allIs')' > minz, :), 50);
                    bar(x, y, 1, 'stack')
            
                    colors = ['brgk'];
                    M = numel(r.dataSet);
                    for n=1:numel(r.dataSet)
                        subplot(4,M,2*M+n)
                        cla
                        scores = r.dataSet(n).meanTypeScores;
                        boundFracs = regionStats(n).czx(:,2)./sum(regionStats(n).czx')';
                        plot(scores, boundFracs, ['o' colors(n)]);
                    end
                    
                    subplot(4,1,4)
                    dispModels(model)
                    drawnow;
                end
            end
            
            % Save to disk every c.trialsPerWrite trials
            if mod(trialNum, r.trialsPerWrite) == 0
                r.Is = Is(1:jotNum,:);
                r.models = models(1:jotNum,:);
                r.numTrials = trialNum;
                r.numJots = jotNum;
                save([c.dir '/' r.fileName], 'r');
            end
        end
    end
    
    % Mark run as finished and exit
    r.Is = Is(1:jotNum,:);
    r.models = models(1:jotNum,:);
    r.numTrials = trialNum;
    r.numJots = jotNum
    save([c.dir '/' r.fileName], 'r');
    
    % Save finished file
    x = 'finished!'
    save([c.dir '/ ' r.finishedFile], 'x');
end
