function d = cohend(matA, matB)
    % Remove NaNs
    matA = matA(~isnan(matA));
    matB = matB(~isnan(matB));
    
    % Sample sizes
    nA = size(matA,1);
    nB = size(matB,1);
    
    % Means
    meanA = mean(matA);
    meanB = mean(matB);
    
    % Sample standard deviations
    stdA = std(matA, 0); % 0 = sample std
    stdB = std(matB, 0);
    
    % Pooled standard deviation
    pooledStd = sqrt(((nA - 1)*stdA^2 + (nB - 1)*stdB^2) / (nA + nB - 2));
    
    % Cohen's d
    d = (meanA - meanB) / pooledStd;
end