function [models, emats] = loadEndModelsAndEmats(directory, N)
    % Load all valid run files
    files = dir(directory);
    j=1;
    disp(['Loading files from ' directory ':']);
    runNums = [];
    for i=1:length(files)
        f = files(i);
        fname = f.name;
        %disp(fname);
        if strfind(fname, 'run_')
            runname = [directory '/' fname];
            try
                load(runname, 'r');
                runs(j) = r;
                runNums(j,1) = r.runNum;
                j=j+1;
            catch
                disp([runname ' is not good.']);
            end
        end
    end
    
    % Sort runs by run number
    A = [runNums, (1:numel(runNums))'];
    B = sortrows(A);
    runs = runs(B(:,2));
    
    % Retrieve last N models from each run
    for i=1:numel(runs)
        models((1+(i-1)*N):(i*N)) = runs(i).models((end-N+1):end);
    end
    
    % Fix each model's gauge
    M = numel(models);
    L = r.L;
    for i=1:M
        models(i) = fixModelGauge(models(i));
    end
    
    % Gather emats together
    emats = zeros(L, 4, M);
    for i=1:M
        emats(:,:,i) = models(i).emat;
    end
end