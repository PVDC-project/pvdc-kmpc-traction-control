function [y,settings] = mapstd_custom(x,varargin)
% custom function for more efficient data scaling
% maps the dataset to zero mean and unit standard deviation
% use cases:
%   I) one input: initial scaling, get mean and std, scale, return settings
%   II) multiple inputs: scale with the provided settings (mean and std)
% the function operates on the rows and matches MATLAB's "mapstd" function
% interface

% I) get mean and std, scale, return settings too
if nargin == 1
    xmean = mean(x,2);
    xstd = std(x,0,2);
    y = (x-xmean) ./ xstd;
    settings.xmean = xmean;
    settings.xstd = xstd;
    return
end

% II) 'apply' or 'reverse' with already calculated settings
if ischar(x) && nargin == 3     % third input are the settings
    operation = x;              % apply or reverse
    x = varargin{1};            % data is the second input
    settings = varargin{2};
    if strcmp(operation,'apply')
        y = (x-settings.xmean) ./ settings.xstd;
    elseif strcmp(operation,'reverse')
        y = x.*settings.xstd + settings.xmean;
    end
    return
end

error('should not have reached this')
end