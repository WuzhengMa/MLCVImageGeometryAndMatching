%Harris interest point detector
imgExample = imread('lenna.pgm');
imshow(checkerboard);
[cim, r, c] = harris(checkerboard, 0.1, 3);

hold on;
plot(r,c, 'yx');
