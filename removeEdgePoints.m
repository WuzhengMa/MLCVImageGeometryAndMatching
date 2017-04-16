function [x, y] = removeEdgePoints(imgExample, x, y, patchSize)
    %Remove edges at boundary
    [R, C] = size(imgExample);
    removalInx = [];
    for i = 1:size(x)
        if(x(i) < patchSize/2 || x(i) > R - patchSize/2)
            removalInx = [removalInx, i];
        end
        if(y(i) < patchSize/2 || y(i) > C - patchSize/2)
            removalInx = [removalInx, i];
        end
    end
    removalInx = unique(removalInx);
    x(removalInx) = [];
    y(removalInx) = [];
end