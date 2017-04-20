init;
patchSize = 32;

%% Q1.2.a) Harris interest point detector
imageName1 = 'FDFigures/DSC02736.ppm';
imageName2 = 'FDFigures/DSC02737.ppm';

if size(size(imread(imageName1)),2) == 2
    imgExample1 = (imread(imageName1));
    imgExample2 = (imread(imageName2));
else 
    imgExample1 = rgb2gray(imread(imageName1));
    imgExample2 = rgb2gray(imread(imageName2));
end


[x1, y1] = harrisDetector(imageName1, patchSize, 1000, 10);
[x2, y2] = harrisDetector(imageName2, patchSize, 1000, 10);


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

%% Q1.3.a) Q1.3.c)
[h, status1] = getHomographyMatrix(x1, y1, x2, y2, 'RANSAC', 2000);

%Obtain the projection points from image 2 to image 1
homoTransPoints = h\[x2;y2;ones(1,size(x2,2))];
oneOverHomoZ=(1./homoTransPoints(3,:));
oneOverHomoZ=[oneOverHomoZ; oneOverHomoZ; oneOverHomoZ];
transPoints = oneOverHomoZ.*homoTransPoints; %Change from homogeneous coordinate to imhomogeneous coordinates

%Calculate the homography accuracy HA
[HA, HD] = getHomoAccuracy([x1; y1], transPoints([1,2], :));

% Draw the projected point back to image 1
imshow(imageName1);
title('Projected figure from image 2 to image 1');
hold on;
plot(transPoints(2,:), transPoints(1,:), 'rx');
hold off;


%% Q1.3.b) Q1.3.d) Computing fundamental matrix F, where x'^TFx = 0
%Using MATLAB build-in function
%[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'Norm8Point');
%[F,inliersIndex] = estimateFundamentalMatrix([x1', y1'],[x2', y2'], 'Method', 'RANSAC')

%Using own implementation
%[F,inliersIndex] = getFundamentalMatrix([x1', y1'],[x2', y2'], 'Norm8Points', 1);
[F,inliersIndex] = getFundamentalMatrix([x1', y1'],[x2', y2'], 'RANSAC', 2000);

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
imtool(stereoAnaglyph(imgExample1,imgExample2));
disparityRange = [0, 96];
disparityMap = disparity(imgExample1, imgExample2, 'BlockSize', 15, 'DisparityRange', disparityRange);
subplot(1,2,2);
imshow(disparityMap,disparityRange);
title('Disparity Map');
colormap jet;
colorbar;

%% 
