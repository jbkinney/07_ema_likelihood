% Quantize data into a data set
function dataSet = quantizeData(data, dataPerBin)
    % Sort scores and indices
    A = [data.regionScores, data.regionIndices];
    A = sortrows(A,-1);
    regionScores = A(:,1);
    regionIndices = A(:,2);
    
    % Quantize
    regionTypes = floor((0:(size(A,1)-1))/dataPerBin)';
    numTypes = numel(unique(regionTypes));
    meanTypeScores = zeros(numTypes,1);
    for i=1:numTypes
        meanTypeScores(i,1) = mean(regionScores(regionTypes == i-1));
    end
    
    % Assign dataSet fields
    dataSet.regionIndices = regionIndices;
    dataSet.regionTypes = regionTypes;
    dataSet.regionScores = regionScores;
    dataSet.numTypes = numTypes;
    dataSet.meanTypeScores = meanTypeScores;
end