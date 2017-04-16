function res = getDescriptors(imgExample, x, y, patchSize, colorHistogram)
    descriptors = zeros(size(x,1), patchSize*patchSize);  % descriptor of size N*M
    colorHistDescriptors = zeros(size(x,1), 256);
    for k = 1:size(x)
        xPatch = x(k) - (patchSize/2);
        yPatch = y(k) - (patchSize/2);
        patchData = zeros(1, patchSize*patchSize);
        for i = 1:patchSize
            for j = 1:patchSize
                patchData(:,(i-1)*patchSize+j) = imgExample(xPatch+i, yPatch+j);
            end
        end
        if colorHistogram
            colorHistDescriptors(k,:) = histcounts(patchData, 256);
        else
            descriptors(k,:) = patchData;
        end
    end
    if colorHistogram
       res = colorHistDescriptors;
    else
       res = descriptors;
    end
end