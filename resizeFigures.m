resize('DSC02715.JPG');
resize('DSC02719.JPG');
resize('DSC02721.JPG');
resize('DSC02722.JPG');
resize('DSC02723.JPG');
resize('DSC02726.JPG');
resize('DSC02728.JPG');
resize('DSC02729.JPG');

function resize(imageName)
    imgExample = imresize(imread(imageName),0.1);
    nameFields = strsplit(imageName, '.');
    name = strcat('HGFigures/', nameFields{1}, '.ppm');
    imwrite(imgExample, name);
end

