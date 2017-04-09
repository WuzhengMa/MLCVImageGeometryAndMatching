%Manually get coordinates of interest points
imgExample = imread('img1.pgm');
imgExample2 = imread('img2.pgm');
%image(imgExample);
imshow(imgExample);
[x,y] = ginput

hold on
for i=1:size(x)
    plot(x(i),y(i), 'yx');
end
