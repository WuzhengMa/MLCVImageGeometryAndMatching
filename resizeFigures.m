resize('IMG_0964.JPG');
resize('IMG_0965.JPG');
resize('IMG_0966.JPG');

function resize(imageName)
    imgExample = imresize(imread(imageName),0.1);
    nameFields = strsplit(imageName, '.');
    name = strcat('FDFigures/', nameFields{1}, '.ppm');
    imwrite(imgExample, name);
end

