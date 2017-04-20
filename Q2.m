init;
patchSize = 32;

%% Q2.1 Find interest points on origin/resized image
imageName1 = 'HGFigures/DSC02723.ppm';
%imageName1ResizedHalfColor = 'HGFigures/DSC02723.ppm';
imageName1ResizedHalfColor = 'HGFigures/DSC02723ResizedHalf.ppm';
%imageName1ResizedHalfColor = (imread(imageName1));
%imageName1ResizedHalfColor = imresize(imageName1ResizedHalfColor, 0.5);
%imwrite(imageName1ResizedHalfColor,'HGFigures/DSC02715ResizedHalf.ppm','ppm')

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


%Show interest points
subplot(2,2,1);
imshow(imageName1);
title('Interest points detection for original image');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName1ResizedHalfColor);
title('Interest points detection for resized image');
hold on;
%plot(features2(1,:), features2(2,:), 'rx');
plot(y2,x2,'rx');
hold off;


%Match interest points
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');

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

%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x2;y2;ones(1,size(x2,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x1(inliers); y1(inliers)], transPoints([1,2], inliers))

%% Q2.1.b) manual v.s. auto
imageName1 = 'HGFigures/DSC02723.ppm';
imageName2 = 'HGFigures/DSC02726.ppm';


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
[features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 5);
[features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3, 'PeakThresh', 5);

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
title('Interest points auto detection for image 1');
hold on;
%plot(features1(1,:), features1(2,:), 'rx');
plot(y1,x1,'rx');
hold off; 

subplot(2,2,2);
imshow(imageName2);
title('Interest auto detection for image 2');
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
[h, inliers, HMstatus] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', 2000);

%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x2;y2;ones(1,size(x2,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x1(inliers); y1(inliers)], transPoints([1,2], inliers))

%Calculate homography from manual selected corrspondence points
[hM, inliersM, HMstatusM] = getHomographyMatrix(x1M, y1M, x2M, y2M, 'RANSAC', 2000);

%Obtain the projection points from image 2 to image 1
homoTransPointsM = hM\[x2M;y2M;ones(1,size(x2M,2))];
oneOverHomoZM=(1./homoTransPointsM(3,:));
oneOverHomoZM=[oneOverHomoZM; oneOverHomoZM; oneOverHomoZM];
transPointsM = oneOverHomoZM.*homoTransPointsM; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HAM, HDM] = getHomoAccuracy([x1M(inliersM); y1M(inliersM)], transPointsM([1,2], inliersM))

%% Q2.1.c) Vary the number of interest points
imageName1 = 'HGFigures/DSC02723.ppm';
imageName2 = 'HGFigures/DSC02726.ppm';


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
[features1, descriptors1] = vl_sift(single(imgExample1), 'EdgeThresh', 3, 'PeakThresh', 5);
[features2, descriptors2] = vl_sift(single(imgExample2), 'EdgeThresh', 3, 'PeakThresh', 5);

%Match interest points by SIFT
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1', descriptors2', features1(1,:)', ...
        features1(2,:)', features2(1,:)', features2(2,:)', 'RANSAC');
    
pointsProvided = ceil(linspace(4, size(x1,2), 10));   % Test by 10 trials
HATrials = zeros(10, 1);
HDTrials = zeros(10, 1);

for i = 1:10    % 10 trials
    %for j = 1:100 % repeat 100 times on every trial
        idxSel = randperm(size(x1,2),pointsProvided(i));
        x1Sel = x1(idxSel);
        y1Sel = y1(idxSel);
        x2Sel = x2(idxSel);
        y2Sel = y2(idxSel);
        %Calculate homography
        [h, inliers, HMstatus] = getHomographyMatrix(x1Sel, y1Sel, x2Sel, y2Sel, 'RANSAC', 2000);

        %Obtain the projection points from image 2 to image 1
        homoTransPoints = h\[x2Sel;y2Sel;ones(1,size(x2Sel,2))];
        oneOverHomoZ=(1./homoTransPoints(3,:));
        oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
        transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

        %Calculate the homography accuracy HA
        [HATrials(i), HDTrials(i)] = getHomoAccuracy([x1Sel; y1Sel], transPoints([1,2], :));

    %end
end

%% Q2.1.c) Vary the number of interest points
imageName1 = 'FDFigures/building1.ppm';
imageName2 = 'FDFigures/building2.ppm';


if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end

%Find interest points of images
[x1, y1] = harrisDetector(imageName1, patchSize, 10000, 10);
[x2, y2] = harrisDetector(imageName2, patchSize, 10000, 10);

%Get color histogram descriptor
colorHistogram = true; %Using the intensity at each pixel of the map is better than using color histogram
descriptors1 = getDescriptors(imgExample1, x1, y1, patchSize, colorHistogram);
descriptors2 = getDescriptors(imgExample2, x2, y2, patchSize, colorHistogram);

%Match interest points by Harris
[descriptors1, descriptors2, x1, y1, x2, y2] = ...
    matchDescriptorSize(descriptors1, descriptors2, x1, y1, x2, y2, 'RANSAC');

%Calculate fundamental matrix
%[F,inliersIndex] = getFundamentalMatrix([x1', y1'],[x2', y2'], 'RANSAC', 20000);
[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'RANSAC', 'NumTrials');

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