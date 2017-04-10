patchSize = 32;

%% Q1.2.a) Harris interest point detector
imageName1 = 'im0.ppm';
imageName2 = 'im1.ppm';

imgExample1 = rgb2gray(imread(imageName1));
[x1, y1] = harrisDetector(imageName1, patchSize);
%coor1 = InterestPointDetector(imgExample1, 0.8, true, 30);
imgExample2 = rgb2gray(imread(imageName2));
%coor2 = InterestPointDetector(imgExample2, 0.8, true, 30);
[x2, y2] = harrisDetector(imageName2, patchSize);



%% Q1.2.b)
colorHistogram = true;
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);
%Show interest points
imshow(imageName1);
title('Test image 1 before rearrange');
hold on;
plot(y1,x1,'rx');
hold off; 

figure;
imshow(imageName2);
title('Test image 2 before rearrange');
hold on;
plot(y2,x2,'rx');
hold off;

[descriptors1, descriptors2, x1, y1, x2, y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2);

%Show interest points
figure;
subplot(2,2,1);
imshow(imageName1);
title('Test image 1 after rearrange');
hold on;
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName2);
title('Test image 2 after rearrange');
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
%[x2, y2] = rearrangePoints(nearestIndex, x2, y2);
[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'Norm8Point');

%% Q1.3.c)
%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x2;y2;ones(1,size(x2,2))];
transPoints = (1./homoTransPoints(3,:)).*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

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
    imgExample = rgb2gray(imread(imageName));
    [cim, x, y] = getAutoInterestPoints(imgExample, 5350, 30);
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

function [x, y] = rearrangePoints(nearestIndex, x, y)
    tempX = x;
    tempY = y;
    %Rearrange x and y to match coordinates in x1 and y1 
    for i = 1:size(nearestIndex)
        x(i) = tempX(nearestIndex(i));
        y(i) = tempY(nearestIndex(i));
    end
end

function [desp1, desp2, X1, Y1, X2, Y2] = matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2)
    %When the size of descriptors are different, keep the pair of descriptors
    %with thee size equals to the descriptor which has smaller size
    %Do nothing when the size of descriptors are the same
    desp1 = [];
    desp2 = [];
    X1 = [];
    Y1 = [];
    X2 = [];
    Y2 = [];
    inlinerFactor = 1.8;
    if size(descriptors2, 1) >= size(descriptors1, 1)
        [nearestIndex, distance] = knnsearch(descriptors2, descriptors1, 'Distance', 'cityblock');
        %Keep only 8 nearest points to reduce error in calculating
        %homography and fundamental matrix
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
                desp1 = [desp1; descriptors1(eightIndex(i),:)];
                desp2 = [desp2; descriptors2(eightNNeighborIndex(i),:)];
                X1 = [X1, x1(eightIndex(i),:)];
                Y1 = [Y1, y1(eightIndex(i),:)];
                X2 = [X2, x2(eightNNeighborIndex(i), :)];
                Y2 = [Y2, y2(eightNNeighborIndex(i), :)];
            end
        end
    elseif size(descriptors2, 1) < size(descriptors1, 1)
        [nearestIndex, distance] = knnsearch(descriptors1, descriptors2, 'Distance', 'cityblock');
        
        %Keep only 8 nearest points to reduce error in calculating
        %homography and fundamental matrix
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
        
        inlinerFlag = distance < mean(distance)*inlinerFactor;
        for i = 1:size(inlinerFlag)
            if inlinerFlag(i)
                desp1 = [desp1; descriptors1(eightNNeighborIndex(i),:)];
                desp2 = [desp2; descriptors2(eightIndex(i),:)];
                X1 = [X1, x1(eightNNeighborIndex(i), :)];
                Y1 = [Y1, y1(eightNNeighborIndex(i), :)];
                X2 = [X2, x2(eightIndex(i),:)];
                Y2 = [Y2, y2(eightIndex(i),:)];
            end
        end
    end
end
