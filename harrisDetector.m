function [x, y] = harrisDetector(imageName, patchSize)
    if size(size(imread(imageName)),2) == 2
        imgExample = (imread(imageName));
    else
        imgExample = rgb2gray(imread(imageName));
    end
    [cim, x, y] = getAutoInterestPoints(imgExample, 2000, 30);
    [x, y] = removeEdgePoints(imgExample, x, y, patchSize);
end
