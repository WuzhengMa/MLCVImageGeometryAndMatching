init;
patchSize = 32;

%% Q1.2.a) Harris interest point detector
imageName1 = 'barn2/im0.ppm';
imageName2 = 'barn2/im1.ppm';

if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

[x1, y1] = harrisDetector(imageName1, patchSize);
[x2, y2] = harrisDetector(imageName2, patchSize);


%% Q1.2.b)
%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

%feature detection and descriptor by SIFT
%[features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3);
%[features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3);

%Show interest points
subplot(1,2,1);
imshow(imageName1);
title('Interest points detection for Test image 1');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(1,2,2);
imshow(imageName2);
title('Interest points detection Test image 2');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2,x2,'rx');
hold off;

%[descriptors1, descriptors2, x1, y1, x2, y2] = matchDescriptorSize(descriptors1', descriptors2', features1(1,:)', features1(2,:)', features2(1,:)', features2(2,:)', 'Norm8Points');
%[descriptors1, descriptors2, x1, y1, x2, y2] = matchDescriptorSize(descriptors1', descriptors2', features1(1,:)', features1(2,:)', features2(1,:)', features2(2,:)', 'RANSAC');
%[descriptors1, descriptors2, x1, y1, x2, y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'Norm8Points');
[descriptors1, descriptors2, x1, y1, x2, y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');

%Show interest points
figure;
subplot(1,2,1);
imshow(imageName1);
title('NN paired patches Test image 1');
hold on;
plot(y1,x1,'rx');
hold off; 

subplot(1,2,2);
imshow(imageName2);
title('NN paired patches Test image 2');
hold on;
plot(y2,x2,'rx');
hold off;


%% Q1.2.c)
%nearestIndex[i] = x means that ith descriptor1 is nearest to xth 
%descriptor2
%[nearestIndex, distance] = knnsearch(descriptors2, descriptors1, 'Distance', 'cityblock');

%% Q1.3.a)
h = getHomographyMatrix(x1, y1, x2, y2);

%% Q1.3.b) Computing fundamental matrix F, where x'^TFx = 0
%[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'Norm8Point');
[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'RANSAC');

%% Q1.3.c)
%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x2;y2;ones(1,size(x2,2))];
transPoints = (1./homoTransPoints(3,:)).*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x1; y1], transPoints([1,2], :))

% Draw the projected point back to image 1
subplot(2,2,[3,4]);
imshow(imageName1);
title('Projected figure from image 2 to image 1');
hold on;
plot(transPoints(2,:), transPoints(1,:), 'rx');
hold off;


%% Q1.3.d)
figure;
subplot(1,2,1);
imshow(imageName1);
title('Epipolar lines of image 1');
hold on;
plot(y1,x1,'go');
%epiLines = epipolarLine(F', [x1, y1]);
epiLines = (F*[x2', y2', ones(size(x2,2),1)]')';
points = lineToBorderPoints(epiLines, [size(imgExample1,2), size(imgExample1,1)]);
line(points(:,[2,4])',points(:,[1,3])');
hold off;

subplot(1,2,2);
imshow(imageName2);
title('Epipolar lines of image 2');
hold on;
plot(y2,x2,'go');
%epiLines = epipolarLine(F, [x2, y2]);
epiLines = (F'*[x1', y1', ones(size(x1,2),1)]')';
points = lineToBorderPoints(epiLines, [size(imgExample1,2), size(imgExample1,1)]);
line(points(:,[2,4])',points(:,[1,3])');
hold off;

%% Utility functions 

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

function [x, y] = harrisDetector(imageName, patchSize)
    if size(size(imread(imageName)),2) == 2
        imgExample = (imread(imageName));
    else
        imgExample = rgb2gray(imread(imageName));
    end
    [cim, x, y] = getAutoInterestPoints(imgExample, 2000, 30);
    [x, y] = removeEdgePoints(imgExample, x, y, patchSize);
end


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

function res = getHomographyMatrix(x1, y1, x2, y2)
    %[x2, y2] = rearrangePoints(nearestIndex, x2, y2);
   
    n = size(x1,2);
%     h = [];
%     for i=1:n
%         rows0 = zeros(3, 1);
%         rowsXY = -[x1(i); y1(i); 1];
%         hx = [rowsXY; rows0; x2(i).*x1(i); x2(i).*y1(i); x2(i)];
%         hy = [rows0; rowsXY; y2(i).*x1(i); y2(i).*y1(i); y2(i)];
%         h = [h, hx, hy]; 
%     end
    rows0 = zeros(3, n);
    rowsXY = -[x1; y1; ones(1,n)];
    hx = [rowsXY; rows0; x2.*x1; x2.*y1; x2];
    hy = [rows0; rowsXY; y2.*x1; y2.*y1; y2];
    h = [hx hy];
    [U, ~, V] = svd(h');
    res = (reshape(V(:,9), 3, 3)).';
end

% function [x, y] = rearrangePoints(nearestIndex, x, y)
%     tempX = x;
%     tempY = y;
%     %Rearrange x and y to match coordinates in x1 and y1 
%     for i = 1:size(nearestIndex)
%         x(i) = tempX(nearestIndex(i));
%         y(i) = tempY(nearestIndex(i));
%     end
% end

function [HA, HD] = getHomoAccuracy(oriPoints, transPoints)
    %HD is computed as the sum of difference between transformed points and
    %the original points
    %HA is computed as the rate where the transformed points is identical
    %to the original points
    HA = 0;
    HD = sum(sum(abs(oriPoints - transPoints)));
    correctMapped = HD < [1;1];
    for i = 1:size(correctMapped,2)
        if correctMapped(1,i) && correctMapped(1,i)
            HA = HA + 1;
        end
    end
    HA = HA / size(correctMapped,2); %Homography Accuracy
end

function [desp1, desp2, X1, Y1, X2, Y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, method)
    %When the method is 'Norm8Points', only 8 pairs of interest points will 
    %be used, hence use NN method to find 8 pairs of closest interest points
    %When the method is 'RANSAC', equal pairs of interest points will be
    %used, hence use NN method to find equal pairs of closest interest
    %points
    inlinerFactor = 1.8;
    if size(descriptors2, 1) >= size(descriptors1, 1)
        if strcmp(method, 'Norm8Points')
            [desp1, desp2, X1, Y1, X2, Y2] = keepEightNearestPoints(descriptors1, descriptors2, x1, y1, x2, y2, inlinerFactor);
        elseif strcmp(method, 'RANSAC')
            [desp1, desp2, X1, Y1, X2, Y2] = keepEqualNumPoints(descriptors1, descriptors2, x1, y1, x2, y2, inlinerFactor);
        end
    elseif size(descriptors2, 1) < size(descriptors1, 1)
        if strcmp(method, 'Norm8Points')
            [desp2, desp1, X2, Y2, X1, Y1] = keepEightNearestPoints(descriptors2, descriptors1, x2, y2, x1, y1, inlinerFactor);
        elseif strcmp(method, 'RANSAC')
            [desp2, desp1, X2, Y2, X1, Y1] = keepEqualNumPoints(descriptors2, descriptors1, x2, y2, x1, y1, inlinerFactor);
        end
    end
end

function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = keepEightNearestPoints(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, inlinerFactor)
    %Use NN search to find 8 pairs of closest interest points
    sSizeDesp = [];
    lSizeDesp = [];
    SX = [];
    SY = [];
    LX = [];
    LY = [];
    [nearestIndex, distance] = knnsearch(lSizeDescriptors, sSizeDescriptors, 'Distance', 'cityblock');
    eightNearestDistances = sort(distance);
    eightNearestDistances = eightNearestDistances(1:8);
    for i = 1:size(eightNearestDistances)
        for j = 1:size(distance)
            if eightNearestDistances(i) == distance(j)
                eightIndex(i) = j;
            end
        end
        eightNNeighborIndex(i) = nearestIndex(eightIndex(i));
    end

    inlinerFlag = eightNearestDistances < mean(eightNearestDistances)*inlinerFactor;
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

function [sSizeDesp, lSizeDesp, SX, SY, LX, LY] = keepEqualNumPoints(sSizeDescriptors, lSizeDescriptors, sx, sy, lx, ly, inlinerFactor)
    %Use NN search to find equal pairs of closest interest points
    sSizeDesp = [];
    lSizeDesp = [];
    SX = [];
    SY = [];
    LX = [];
    LY = [];
    [nearestIndex, distance] = knnsearch(lSizeDescriptors, sSizeDescriptors, 'Distance', 'cityblock');
    inlinerFlag = distance < mean(distance)*inlinerFactor;
    for i = 1:size(inlinerFlag)
        if inlinerFlag(i)
            sSizeDesp = sSizeDescriptors;
            lSizeDesp = [lSizeDesp; lSizeDescriptors(nearestIndex(i),:)];
            SX = [SX, sx(i,:)'];
            SY = [SY, sy(i,:)'];
            LX = [LX, lx(nearestIndex(i), :)];
            LY = [LY, ly(nearestIndex(i), :)];
        end
    end
end
