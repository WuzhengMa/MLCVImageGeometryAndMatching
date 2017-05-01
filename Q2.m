init;
patchSize = 32;

%% Q2.1 Find interest points on origin/resized image
imageName1 = 'HGFigures/trashChairs1.ppm';
%imageName1ResizedHalfColor = 'HGFigures/DSC02723.ppm';
imageName1ResizedHalfColor = 'HGFigures/trashChairs1ResizedHalf.ppm';
%imageName1ResizedHalfColor = (imread(imageName1));
%imageName1ResizedHalfColor = imresize(imageName1ResizedHalfColor, 0.5);
%imwrite(imageName1ResizedHalfColor,'HGFigures/trashChairs1ResizedHalf.ppm','ppm')

if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample1ResizedHalf = (imread(imageName1ResizedHalfColor));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample1ResizedHalf = rgb2gray(imread(imageName1ResizedHalfColor));
end

%imshow(imgExample1);

%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize, 500, 10);
[x2, y2] = harrisDetector(imageName1ResizedHalfColor, patchSize, 500, 10);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample1ResizedHalf, x2, y2, patchSize, colorHistogram);

%feature detection and descriptor by SIFT
[features1, descriptors1SIFT] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 5);
[features2, descriptors2SIFT] = vl_sift(single(imgExample1ResizedHalf), 'EdgeThresh', 3, 'PeakThresh', 5);
x1SIFT = features1(1,:);
y1SIFT = features1(2,:);
x2SIFT = features2(1,:);
y2SIFT = features2(2,:);

%Show interest points
subplot(2,2,1);
imshow(imageName1);
title('Interest points detection for original image w/ Harris');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image w/ Harris');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2,x2,'rx');
hold off;
%{
%Show interest points
subplot(2,2,3);
imshow(imageName1);
title('Interest points detection for original image w/ SIFT');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
%plot(y1SIFT,x1SIFT,'rx');
vl_plotframe(features1);
%vl_plotsiftdescriptor(descriptors1SIFT, features1);
hold off; 

subplot(2,2,4);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image w/ SIFT');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
%plot(y2SIFT,x2SIFT,'rx');
vl_plotframe(features2);
%vl_plotsiftdescriptor(descriptors2SIFT, features2);
hold off;
%}
%Match interest points
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');

%[~, ~, x1, y1, x2, y2] = matchSIFTDescriptors(descriptors1SIFT, descriptors2SIFT,...
%            features1(1,:), features1(2,:), features2(1,:), features2(2,:)); 

%Show interest points
subplot(2,2,3);
imshow(imageName1);
title('Interest points detection for original image');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,4);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2,x2,'rx');
hold off;

%Calculate homography
[h, inliers, HMstatus] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', 2000);

%{
%Obtain the projection points from image 2 to image 1
homoTransPoints = h * [x1;y1;ones(1,size(x1,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x2(inliers); y2(inliers)], transPoints([1,2], inliers))
%}
%Calculate HA error
load('Results/GT_Q2_1_a.mat');
[HA, HD] = getHAFromGT(h, GTx1, GTy1, GTx2, GTy2)

%% Q2.1.b) manual v.s. auto
imageName1 = 'HGFigures/trashChairs1.ppm';
imageName2 = 'HGFigures/trashChairs2.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

%{
%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize, 500, 10);
[x2, y2] = harrisDetector(imageName2, patchSize, 500, 10);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

%Match interest points by Harris
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');
%}


%feature detection and descriptor by SIFT
[features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 10);
[features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3, 'PeakThresh', 10);

%Match interest points by SIFT
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1', descriptors2', features1(1,:)', ...
        features1(2,:)', features2(1,:)', features2(2,:)', 'RANSAC');


[x1M, y1M] = getManualInterestPoints(imageName1);
[x2M, y2M] = getManualInterestPoints(imageName2);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1M = getDescriptors(imgExample1, x1M, y1M, patchSize, colorHistogram);
descriptors2M = getDescriptors(imgExample2, x2M, y2M, patchSize, colorHistogram);

%Match interest points
[descriptors1M, descriptors2M, x1M, y1M, x2M, y2M] = ...
    matchDescriptorSize(descriptors1M, descriptors2M, x1M, y1M, x2M, y2M, 'RANSAC');

%Show interest points
subplot(2,2,1);
imshow(imageName1);
title('Interest points by SIFT for image 1');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName2);
title('Interest points by SIFT for image 2');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2,x2,'rx');
hold off;

subplot(2,2,3);
imshow(imageName1);
title('Interest points manually selected for image 1');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1M,x1M,'rx');
hold off; 

subplot(2,2,4);
imshow(imageName2);
title('Interest points manually selected for image 2');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2M,x2M,'rx');
hold off;

%Calculate homography
[h, inliers, HMstatus] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', 20000);

%{
%Obtain the projection points from image 2 to image 1
homoTransPoints = h * [x1;y1;ones(1,size(x1,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x2(inliers); y2(inliers)], transPoints([1,2], inliers))
%}
%load('Results/Q2_1_b.mat');
[HA, HD] = getHAFromGT(h, GTx1, GTy1, GTx2, GTy2)

%Calculate homography from manual selected corrspondence points
[hM, inliersM, HMstatusM] = getHomographyMatrix(x1M, y1M, x2M, y2M, 'RANSAC', 2000);

%{
%Obtain the projection points from image 2 to image 1
homoTransPointsM = hM * [x1M;y1M;ones(1,size(x1M,2))];
oneOverHomoZM=(1./homoTransPointsM(3,:));
oneOverHomoZM=[oneOverHomoZM; oneOverHomoZM; oneOverHomoZM];
transPointsM = oneOverHomoZM.*homoTransPointsM; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HAM, HDM] = getHomoAccuracy([x2M(inliersM); y2M(inliersM)], transPointsM([1,2], inliersM))
%}
[HAM, HDM] = getHAFromGT(hM, GTx1, GTy1, GTx2, GTy2)

%% Q2.1.c) Vary the number of interest points
imageName1 = 'HGFigures/trashChairs1.ppm';
imageName2 = 'HGFigures/trashChairs2.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

%{
%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize, 500, 10);
[x2, y2] = harrisDetector(imageName2, patchSize, 500, 10);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

%Match interest points by Harris
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');
%}


%feature detection and descriptor by SIFT
[features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 10);
[features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3, 'PeakThresh', 10);

%Match interest points by SIFT
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1', descriptors2', features1(1,:)', ...
        features1(2,:)', features2(1,:)', features2(2,:)', 'RANSAC');
    
pointsProvided = ceil(linspace(4, size(x1,2), 10));   % Test by 10 trials
HATrials = zeros(10, 100);
HDTrials = zeros(10, 100);

for i = 1:10    % 10 trials
    for j = 1:100 % repeat 100 times on every trial
        idxSel = randperm(size(x1,2),pointsProvided(i));
        x1Sel = x1(idxSel);
        y1Sel = y1(idxSel);
        x2Sel = x2(idxSel);
        y2Sel = y2(idxSel);
        %Calculate homography
        [h, inliers, HMstatus] = getHomographyMatrix(x1Sel, y1Sel, x2Sel, y2Sel, 'RANSAC', 2000);

        %{
        %Obtain the projection points from image 2 to image 1
        homoTransPoints = h * [x1Sel;y1Sel;ones(1,size(x1Sel,2))];
        oneOverHomoZ=(1./homoTransPoints(3,:));
        oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
        transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

        %Calculate the homography accuracy HA
        [HATrials(i), HDTrials(i)] = getHomoAccuracy([x2Sel; y2Sel], transPoints([1,2], :));
        %}
        
        load('Results/GT_Q2_1_c.mat');
        [HATrials(i,j), HDTrials(i,j)] = getHAFromGT(h, GTx1, GTy1, GTx2, GTy2);
        
    end
end

%% Q2.1.c) Vary the number of interest points
imageName1 = 'FDFigures/trashChairBoard1.ppm';
imageName2 = 'FDFigures/trashChairBoard2.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize, 10, 10);
[x2, y2] = harrisDetector(imageName2, patchSize, 10, 10);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

%Match interest points by Harris
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');

%{
%Calculate homography
[h, inliers, HMstatus] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', 2000);

%Obtain the projection points from image 2 to image 1
q = h * [x1; y1; ones(1, size(x1,2))];
p = q(3,:);
transPoints = [q(1,:)./p; q(2,:)./p];

%Show interest points
subplot(1,2,1);
imshow(imageName1);
title('Interest points auto detection for image 1');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(1,2,2);
imshow(imageName2);
title('Interest auto detection for image 2');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(transPoints(2,:),transPoints(1,:),'rx');
hold off;
%}
    
%Calculate fundamental matrix
[F,inliersIndex] = getFundamentalMatrix([x1', y1'],[x2', y2'], 'RANSAC', 20000);
%[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], ...
%    'Method', 'RANSAC', 'NumTrials', 10000, 'DistanceThreshold', 10);

%Get epipoles
[e1, e2] = getEpipoles(F);

figure;
subplot(1,2,1);
imshow(imageName1);
title('Epipolar lines of image 1');
hold on;
plot(y1,x1,'go');
%epiLines = epipolarLine(F, [x2, y2]);
epiLines = (F*[x2', y2', ones(size(x2,2),1)]')';
points = lineToBorderPoints(epiLines, [size(imgExample1,2), size(imgExample1,1)]);
line(points(:,[2,4])',points(:,[1,3])');
hold off;

subplot(1,2,2);
imshow(imageName2);
title('Epipolar lines of image 2');
hold on;
plot(y2,x2,'go');
%epiLines = epipolarLine(F', [x1, y1]');
epiLines = (F'*[x1', y1', ones(size(x1,2),1)]')';
points = lineToBorderPoints(epiLines, [size(imgExample1,2), size(imgExample1,1)]);
line(points(:,[2,4])',points(:,[1,3])');
hold off;

%% Q2.2.C) Compute the disparity map                             
figure;
subplot(1,2,1);
title('Red-cyan graph of two images');
imshow(stereoAnaglyph(imgExample1,imgExample2));
%imtool(stereoAnaglyph(imgExample1,imgExample2));
disparityRange = [10, 58]; %[1, 97]
disparityMap = disparity(imgExample1, imgExample2, 'BlockSize', 15, 'DisparityRange', disparityRange);
subplot(1,2,2);
imshow(disparityMap,disparityRange);
title('Disparity Map');
colormap jet;
colorbar;

%% Q2.2.D) Compute the depth map
%{
focalLen = 0.018;   %18mm
baseline = 0.1;     %10cm
depthMap = (focalLen * baseline) ./ disparityMap; 
[gridX, gridY] = meshgrid(1:size(disparityMap,2), 1:size(disparityMap, 1));
ptCloud = pointCloud(cat(3, gridX, gridY, disparityMap), 'Color', (imread(imageName1)));

% Create a streaming point cloud viewer
player3D = pcplayer([1, 400], [1, 300], [0, 100], 'VerticalAxis', 'y', ...
    'VerticalAxisDir', 'down');

% Visualize the point cloud
view(player3D, ptCloud);
%}

focalLen = 18;   %18mm
baseline = 100;     %10cm
depthMap = (focalLen * baseline) ./ disparityMap; 
[gridX, gridY] = meshgrid(1:size(disparityMap,2), 1:size(disparityMap, 1));
ptCloud = pointCloud(cat(3, gridX, gridY, depthMap), 'Color', (imread(imageName1)));

% Create a streaming point cloud viewer
player3D = pcplayer([1, 400], [1, 300], [min(min(depthMap)), max(max(depthMap))], 'VerticalAxis', 'y', ...
    'VerticalAxisDir', 'down');


% Visualize the point cloud
view(player3D, ptCloud);


%% Q2.2.E) Change the focal length and compute the depth map
imageName1 = 'FDFigures/trashChairBoard1_zoomin.ppm';
imageName2 = 'FDFigures/trashChairBoard2_zoomin.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

disparityRange = [10, 58];
disparityMap = disparity(imgExample1, imgExample2, 'BlockSize', 15, 'DisparityRange', disparityRange);

%{
focalLen = 0.020;   %20mm
baseline = 0.1;     %10cm
depthMap = (focalLen * baseline) ./ disparityMap; 
[gridX, gridY] = meshgrid(1:size(disparityMap,2), 1:size(disparityMap, 1));
ptCloud = pointCloud(cat(3, gridX, gridY, disparityMap), 'Color', (imread(imageName1)));
%}

focalLen = 18.5;   %18.5mm
baseline = 100;     %10cm
depthMap = (focalLen * baseline) ./ disparityMap; 
[gridX, gridY] = meshgrid(1:size(disparityMap,2), 1:size(disparityMap, 1));
ptCloud = pointCloud(cat(3, gridX, gridY, depthMap), 'Color', (imread(imageName1)));

% Create a streaming point cloud viewer
player3D = pcplayer([1, 400], [1, 300], [min(min(depthMap)), max(max(depthMap))], 'VerticalAxis', 'y', ...
    'VerticalAxisDir', 'down');

% Visualize the point cloud
view(player3D, ptCloud);

%% Q2.2.E) Add noise to disparity and compute the depth map

imageName1 = 'FDFigures/trashChairBoard1.ppm';
imageName2 = 'FDFigures/trashChairBoard2.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

disparityRange = [10, 58];
disparityMap = disparity(imgExample1, imgExample2, 'BlockSize', 15, 'DisparityRange', disparityRange);


disparityMapRng = disparityMap + normrnd(0, 2/3, [size(disparityMap)]);
%disparityMapRng = disparityMapRng ./ (max(max(disparityMapRng))-min(min(disparityMapRng)))...
%    .* (disparityRange(2)-disparityRange(1));
%disparityMapRng = disparityRange(1)-min(min(disparityMapRng)) + disparityMapRng;

focalLen = 18;   %18mm
baseline = 100;     %10cm
depthMap = (focalLen * baseline) ./ disparityMapRng; 
[gridX, gridY] = meshgrid(1:size(disparityMapRng,2), 1:size(disparityMapRng, 1));
ptCloud = pointCloud(cat(3, gridX, gridY, depthMap), 'Color', (imread(imageName1)));

% Create a streaming point cloud viewer
player3D = pcplayer([1, 400], [1, 300], [min(min(depthMap)), max(max(depthMap))], 'VerticalAxis', 'y', ...
    'VerticalAxisDir', 'down');

% Visualize the point cloud
view(player3D, ptCloud);