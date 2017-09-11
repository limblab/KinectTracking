function [ out ] = inBound( val, window )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    out = val >window(1) & val <window(2);
end

