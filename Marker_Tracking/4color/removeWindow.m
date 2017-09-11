function [ color2 ] = removeWindow( color1, indices, xWin, yWin, zWin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    color2 = color1;
    for i = 1:length(indices)
        temp=color1{indices(i)};
        x1=temp(1:end/3);
        y1=temp(end/3+1:2*end/3);
        z1=temp(2*end/3+1:end);
        pixels = [x1;y1;z1]';
        pixels((inBound(x1, xWin) & inBound(y1, yWin) & inBound(z1,zWin)),:) = []; 
        color2{indices(i)} = [pixels(:,1); pixels(:,2); pixels(:,3)]';
    end
end

