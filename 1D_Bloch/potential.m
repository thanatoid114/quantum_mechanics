function [ output_args ] = potential( x,v1,a1,a2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
output_args = zeros(size(x));

output_args = output_args+((x>=a1)&(x<= a2))*v1;
    

end

