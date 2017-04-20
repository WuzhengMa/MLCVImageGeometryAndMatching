function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = keepEqualNumPointsSortedDist(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, inlinerFactor)
    %Use NN search to find equal pairs of closest interest points
    sSizeDesp = [];
    lSizeDesp = [];
    SX = [];
    SY = [];
    LX = [];
    LY = [];
    dist = [];
    [nearestIndex, distance] = knnsearch(lSizeDescriptors, sSizeDescriptors, 'Distance', 'cityblock');
    inlinerFlag = distance < mean(distance)*inlinerFactor;
    for i = 1:size(inlinerFlag)
        if inlinerFlag(i)
            sSizeDesp = [sSizeDesp; sSizeDescriptors(i,:)];
            lSizeDesp = [lSizeDesp; lSizeDescriptors(nearestIndex(i),:)];
            SX = [SX, sx(i,:)'];
            SY = [SY, sy(i,:)'];
            LX = [LX, lx(nearestIndex(i), :)];
            LY = [LY, ly(nearestIndex(i), :)];
            dist = [dist distance(i)];
        end
    end

    [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = sortDistance(sSizeDesp, lSizeDesp, SX, SY, LX, LY, dist);

function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = sortDistance(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, dist)
    [~, idx]=sort(dist);
    M = size(sSizeDescriptors, 2); %size of a single desc
    N = size(sSizeDescriptors, 1); %No of descps
    sSizeDesp = zeros(N, M);
    lSizeDesp = zeros(N, M);
    SX = zeros(1, N);
    SY = zeros(1, N);
    LX = zeros(1, N);
    LY = zeros(1, N);
    for i = 1:N
        sSizeDesp(i,:) = sSizeDescriptors(idx(i),:);
        lSizeDesp(i,:) = lSizeDescriptors(idx(i),:);
        SX(:,i) = sx(idx(i));
        SY(:,i) = sy(idx(i));
        LX(:,i) = lx(idx(i));
        LY(:,i) = ly(idx(i));
    end

    