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