% Madeline Sprague (2/1/2024) - mjsprague@uri.edu
%
% This script models a source-receiver scenario using inputted information
% and HYCOM output to find the expected arrival time of a signal produced
% by the source (position 1) and received at position 2. 
%
% Steps: 
%   user inputs lat/lon coordinates of positions 1 and 2, along with
%   source/receiver depths. hycom_bh models the propagation, specifying a 
%   receiver at position 2. The script will determine the 'main' arrival at 
%   the receiver by the strongest cluster of arrival times. The script 
%   outputs an arrival time with an associated uncertainty. 
%
% INPUTS: 
%   lon1, lat1, lon2, lat2: double
%       Decimal lat/lon values associated with source (lat1, lon1) and
%       receiver (lat2, lon2)
%   time: datenum
%       The time that will be used to generate the HYCOM model output. In
%       addition to a datenum value, one can instead use 'summer' or
%       'winter' for the climatology files (only if you have these). 
%   sz, rz: double
%       Source depth and receiver depth in meters. 
%
% OUTPUTS: 
%   atime: double
%       The arrival time, in seconds, predicted by the propagation. 
%   atime_d: double
%       The uncertainty of the arrival time (seconds), solved by the
%       standard deviation of arrival time from the points that make up the
%       main cluster. 
  

function [atime, atime_d] = expected_atime(lon1, lat1, lon2, lat2, time, sz, rz)
    
% generate struct to input user settings without manual inputs 

    rx = m_lldist([lon1 lon2], [lat1 lat2]); 
    outputs = struct('filename',   'atime_temp', ...
                     'env_title',  'atime_temp', ...
                     'runtype',    '''A''',      ... % arrivals mode 
                     'freq',       262,          ... % assume f = 262 Hz
                     'dep_source', sz,           ...
                     'nrays',      0,            ... 
                     'fr_min',     -88,          ...
                     'fr_max',     88,           ...
                     'n_rec_z',    1,            ...
                     'z_rec_min',  rz,           ...
                     'z_rec_max',  rz,           ...
                     'n_rec_x',    1,            ...
                     'x_rec_min',  rx,           ...
                     'x_rec_max',  rx,           ...
                     'cp_bot',     1800,         ...
                     'cs_bot',     400,          ...
                     'rho_bot',    1800,         ...
                     'apt',        0.104,        ...
                     'ast',        0.26); 
                 
% we model the receiver at lon2, lat2, which shouldn't be at the very edge
% of the propagation - so extend the lat/lon transect just a little 

    dist_lat = lat2 - lat1; 
    dist_lon = lon2 - lon1; 
    lat2     = lat1 + dist_lat * 1.05; 
    lon2     = lon1 + dist_lon * 1.05; 
                     
% run hycom_bh 

    fpath = cd; 
    [outputs] = hycom_bh(lat1, lon1, lat2, lon2, time, outputs); 
    close all 

% analyze arrivals data 

    arrfile = [fpath '\' outputs.filename '\' outputs.filename '.arr']; 
    [arr, ~] = read_arrivals_asc(arrfile); 
    atimes = double(arr.delay); 
    atimes_hist = floor(min(atimes)):ceil(max(atimes)); 
    atimes_density = NaN(size(atimes_hist)); 
    for i = 1:length(atimes_hist)-1
        atimes_density(i) = sum(atimes > atimes_hist(i) & atimes < atimes_hist(i+1)); 
    end
    atime_rough = atimes_hist(atimes_density == max(atimes_density)); 
        atime_rough = mean(atime_rough); 
    idx         = atimes > atime_rough - 0.5 & atimes < atime_rough + 1.5; 
    atime       = mean(atimes(idx)); % final value (s)
    atime_d     = std(atimes(idx));  % associated uncertainty (s)

% delete files 

    rmdir(outputs.filename, 's'); 

end