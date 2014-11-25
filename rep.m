function [result]=rep(array, count) 
matrix = repmat(array, count,1);
result = matrix(:);

%see link for description: http://stackoverflow.com/questions/14615305/a-similar-function-to-rs-rep-i
%n-matlab