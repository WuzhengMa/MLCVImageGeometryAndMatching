resize('DSC02704.JPG');
resize('DSC02706.JPG');
resize('DSC02707.JPG');
resize('DSC02708.JPG');
resize('DSC02709.JPG');
resize('DSC02710.JPG');
resize('DSC02711.JPG');
resize('DSC02712.JPG');
resize('DSC02713.JPG');
resize('DSC02714.JPG');

function resize(imageName)
    imgExample = imresize(imread(imageName),0.1);
    nameFields = strsplit(imageName, '.');
    name = strcat('FDFigures/', nameFields{1}, '.ppm');
    imwrite(imgExample, name);
end

