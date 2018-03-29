function [x,y,z] = getXYZfromKinect(colorFrame)
    x=colorFrame(1:end/3);
    y=colorFrame(end/3+1:2*end/3);
    z=colorFrame(2*end/3+1:end);
