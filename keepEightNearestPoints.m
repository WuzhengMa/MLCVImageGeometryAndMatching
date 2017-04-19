function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = keepEightNearestPoints(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, inlinerFactor)
    %Use NN search to find 8 pairs of closest interest points
    sSizeDesp = [];
    lSizeDesp = [];
    SX = [];
    SY = [];
    LX = [];
    LY = [];
    [nearestIndex, distance] = knnsearch(lSizeDescriptors, sSizeDescriptors, 'Distance', 'cityblock');
    [eightNearestDistances, sortedIdx] = sort(distance);
    eightNearestDistances = eightNearestDistances(1:8);
    eightIndex = sortedIdx(1:8);
    for i = 1:size(eightNearestDistances)
        eightNNeighborIndex(i) = nearestIndex(eightIndex(i));
    end
    %{
    for i = 1:size(eightNearestDistances)
        for j = 1:size(distance)
            if eightNearestDistances(i) == distance(j)
                eightIndex(i) = j;
            end
        end
        eightNNeighborIndex(i) = nearestIndex(eightIndex(i));
    end
    %}

    inlinerFlag = eightNearestDistances <= mean(eightNearestDistances)*inlinerFactor;
    for i = 1:size(inlinerFlag)
        if inlinerFlag(i)
            sSizeDesp = [sSizeDesp; sSizeDescriptors(eightIndex(i),:)];
            lSizeDesp = [lSizeDesp; lSizeDescriptors(eightNNeighborIndex(i),:)];
            SX = [SX, sx(eightIndex(i),:)];
            SY = [SY, sy(eightIndex(i),:)];
            LX = [LX, lx(eightNNeighborIndex(i), :)];
            LY = [LY, ly(eightNNeighborIndex(i), :)];
        end
    end
end