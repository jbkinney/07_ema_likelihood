function quecomputations(computations, runlocally)
delay = 60*60;   % Copy every 15 minutes

% Load computation and set running parameters
for i=1:numel(computations)
    computations(i).runlocally = runlocally;
    computations(i).finished = false;
end

% Run computations serially
for compNum = 1:length(computations)
    
    c = computations(compNum)
    system(['rm -rf ' c.dir]);
    system(['mkdir ' c.dir]);
    save([c.dir '/computation.mat'], 'c');
    
    % If running locally on laci, do only first mcmc run
    if runlocally
        c.main(1, runlocally);       
        
    % If running remotely on fafner , submit all jobs
    else
        
        % Create string containing job-executing script
        jobString = [...
            'runNum = str2num(getenv(''SGE_TASK_ID'')); '...
            'rand(''seed'', runNum); '...
            'runDir = getenv(''RUN_DIR''); '...
            'load([runDir ''/computation.mat'']); '...
            'c.main(runNum, false); '...
            'exit;'];
        
        % Create job command
        f = fopen('job.sh', 'w');
        fprintf(f, ['matlab -nodesktop -nojvm -nosplash -nodisplay -r "' jobString '"']);
        fclose(f);
        system('chmod +x job.sh');

        % Submit jobs (now or not at all)
        % Have send email updates profusely
        system(['qsub' ...
            ' -l 1day' ...
            ' -sync n' ...
            ' -t 1-' num2str(c.numRuns) ':1' ...
            ' -M ' c.recipients ...
            ' -m beas' ...
            ' -v RUN_DIR=''' c.dir '''' ...
            ' -N ''' c.description '''' ...
            ' -cwd -now n ./job.sh']);
        disp([num2str(c.numRuns) ' jobs submitted!']);
        
        % Mark computation as not done. 
        computations(compNum).finished = false;
    end 
end

% Loop as long as any computation is not finished (if not running locally)
allDone = false;
while (~allDone) & (~runlocally)
        
    % Iterate through all computations
    pause(delay);
    allDone = true;
    for compNum = 1:length(computations)
        c = computations(compNum);

        % Copy over results if not finished
        if ~c.finished
            allDone = false;
            
            % Copy contents over
            disp(['Copying contents of ' c.dir ' to laci.'])
            system(['scp -r ' c.dir ' ' c.laciDir '/.']);
            
            % Test to see if finished now
            done = true;
            for n=1:length(c.numRuns)
                done = done & (exist([c.dir '/' c.runs(n).finishedFile], 'file') == 2);
            end
                
            % Record whether or not computation is now finished
            computations(compNum).finished = done;
        end
    end
end

if ~runlocally
    % Copy data one last time
    pause(delay)
    for compNum = 1:length(computations)
        c = computations(compNum);
        disp(['Copying contents of ' c.dir ' to laci.'])
        system(['scp -r ' c.dir ' ' c.laciDir '/.']);
    end
end
