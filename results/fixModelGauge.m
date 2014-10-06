function model = fixModelGauge(model)
    % Bring model into gauge. Determine consensus site
    % and set maximum cutoff to 1. Make sure all cutoffs are > 0.
    % Force all elements <= 1;
    L = size(model.emat,1);
    for i=1:L
        [x,b] = min(model.emat(i,:)');
        model.cutoff = model.cutoff-x;
        model.emat(i,:) = model.emat(i,:) - x;
        model.consensus(i) = b;
    end
    [maxCutoff, model.fixedCutoffNum] = max(model.cutoff(:));
    model.emat = model.emat/maxCutoff;
    model.cutoff = model.cutoff/maxCutoff;
    for i=1:L
        for b=1:4
            model.emat(i,b) = min(1, model.emat(i,b));
        end
    end
    if ~all(model.cutoff(:) > 0)
        error('ERROR: gauge-fixed cutoffs not all > 0')
    end
end