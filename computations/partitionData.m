% Partitions data into two disjoint data sets
function [dataSetA, dataSetB] = partitionData(data, dataPerBin)
    N = numel(data.regionIndices);
    indices = randsample(N,N);
    indicesA = indices(1:round(N/2));
    indicesB = indices((round(N/2)+1):end);
    
    dataA.regionIndices = data.regionIndices(indicesA);
    dataA.regionScores = data.regionScores(indicesA);
    dataSetA = quantizeData(dataA, dataPerBin);
    
    dataB.regionIndices = data.regionIndices(indicesB);
    dataB.regionScores = data.regionScores(indicesB);
    dataSetB = quantizeData(dataB, dataPerBin);
end
    
    
    
    