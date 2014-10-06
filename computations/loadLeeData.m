function [data, A] = loadChIPData(tfName)
    
    % Open file
    seqf = fopen('data/leeData.txt');
    
    % Get number of transcription factors
    C = textscan(seqf, '%s %n', 1, 'delimiter', '\t', 'bufSize', 10000);
    ntfs = C{2};
    
    % Get list of all tf names
    C = textscan(seqf, repmat(['%s '], 1, ntfs+1), 1, 'delimiter', '\t', 'bufSize', 1000000);
    tfNames = {};
    for i=2:length(C)
        tfNames{i-1} = upper(C{i}{1});
    end
    
    % Get list of all TF pvalues
    formatString = ['%s ' repmat(['%n '],1,ntfs) '%s'];
    C = textscan(seqf, formatString, 'delimiter', '\t', 'bufSize', 1000000);
    
    % If TF name not provided, return list of possible TFs
    if ~tfName
        dataSet = [];
        data = tfNames;
        regionIndices = [];
        regionTypes = [];
        for tfIndex=1:ntfs
            pvalues = C{tfIndex+1};
            medpvalue = nanmedian(pvalues);
            numPos(tfIndex,1) = sum(pvalues < pvalueThreshold);
            numNeg(tfIndex,1) = sum(pvalues > medpvalue);
        end
        return
    else
        tfName = upper(tfName);
    end
    
    % Find desired tf
    tfIndex = strmatch(tfName, tfNames, 'exact'); % require exact match
    if ~tfIndex
        error(['ERROR! Unable to find transcription factor ' tfName]);
    end
    
    % Load region names, pvalues and sequences
    names = C{1};
    pvalues = C{tfIndex+1};
    seqs = C{end};
    
    % Get median pvalue
    medpvalue = nanmedian(pvalues(:));
    
    % Store in structure of regions
    nregions = length(names);
    regions = repmat(struct(), nregions, 1);
    regionIndices = NaN(nregions,1);
    for i=1:nregions
        regions(i).tf = tfName;
        regions(i).name = names{i};
        regions(i).pvalue = pvalues(i);
        seq = seqs{i};
        regions(i).seq = seq;
        regions(i).nseq = double(1*(seq=='A') + 2*(seq=='C') + 3*(seq=='G') + 4*(seq=='T'));
    end 
    
    A = [-norminv(pvalues,0,1), (1:nregions)', pvalues];   % Concatenate lrs with indices
    A = A(~isnan(A(:,1)),:);      % Remove rows with NaN lrs
    A = sortrows(A,-1);
    data.regionScores = A(:,1);
    data.regionIndices = A(:,2);
    data.regions = regions;
end