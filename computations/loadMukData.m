% Loads regions from file
function [data, A] = loadMukData(tf)
    fname = 'data/mukData.txt';
    seqf = fopen(fname);
    
	% Read off header line
    D = textscan(seqf, '%s %s %s %s %s %s %s %s %s', 1, 'bufSize', 10000);

	% Read actual data
	C = textscan(seqf, '%s %d %n %n %n %n %n %n %s', 'bufSize', 1000000);
    names = C{1};
    
    tf = lower(tf);
    
    % Return structure array of regions
    nregions = length(C{1});
    regions = repmat(struct(), nregions, 1);
    regionIndices = zeros(nregions,1);
    lrs = zeros(nregions,1);
    for i=1:nregions
        
        regions(i).name = C{1}{i};
        regions(i).length = C{2}(i);
        regions(i).MIG1lr = C{3}(i);
        regions(i).MIG1p = C{4}(i);
        regions(i).ABF1lr = C{5}(i);
        regions(i).ABF1p = C{6}(i);
        regions(i).RAP1lr = C{7}(i);
        regions(i).RAP1p = C{8}(i);
        
        seq = upper(C{9}{i});
        regions(i).seq = seq;
        regions(i).nseq = double(1*(seq=='A') + 2*(seq=='C') + 3*(seq=='G') + 4*(seq=='T'));
    
        if strcmp(tf, 'abf1')
            regions(i).pvalue = regions(i).ABF1p;
            regions(i).lr = regions(i).ABF1lr;
        elseif strcmp(tf, 'rap1')
            regions(i).pvalue = regions(i).RAP1p;
            regions(i).lr = regions(i).RAP1lr;
        elseif strcmp(tf, 'mig1')
            regions(i).pvalue = regions(i).MIG1p;
            regions(i).lr = regions(i).MIG1lr;
        else
            error('Invalid tf');
        end
        
        lrs(i,1) = regions(i).lr;
    end

    A = [lrs, (1:nregions)'];       % Concatenate lrs with indices
    A = A(~isnan(A(:,1)),:);        % Remove rows with NaN lrs
    A = sortrows(A,-1);
    data.regionScores = A(:,1);
    data.regionIndices = A(:,2);
    data.regions = regions;
end