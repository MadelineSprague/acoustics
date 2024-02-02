% Madeline Sprague (2/1/2024) - mjsprague@uri.edu
%
% [bathy, xfill, yfill] = gebco_bathy(lat, lon, x) gets a slice of 
% bathymetry for a lat/lon transect using the GEBCO global elevation grid.
%
% This function makes use of the local file GEBCO_2021_sub_ice_topo.nc,
% which is included with my functions at
% https://github.com/MadelineSprague/acoustics. 
%
% This produces subsurface bathymetry with depth positive downwards, which
% can be edited on line 63. 
% 
% Inputs: 
%   lat, lon: double array 
%       Vectors of lat and lon values for the desired transect.
%   x: double array 
%       The x-axis variable plotted on the transect (e.g. cumulative
%       distance traveled, or it could be lat or lon). Must be a vector of
%       the same size as lat and lon. 
% 
% Outputs: 
%   bathy: double array 
%       Vector of bathymetry depths for the transect, positive downwards.
%   xfill, yfill: double array 
%       Two vectors outlining the fill area if you want to use the fill
%       function to add the bathymetry to your plot. 
%
% (Yes, you may not need xfill or yfill, in which case you can set these
% outputs to ~ and enter [] for x.)

function [bathy, xfill, yfill] = gebco_bathy(lat, lon, x)

% check dimensions

    if size(lat) ~= size(lon) | size(lon) ~= size(x)
        error('Must specify all variables as equal-size vectors')
    end
    
    if isrow(lat)
        lat = lat'; 
        lon = lon'; 
        x   = x'; 
    end

fname     = 'GEBCO_2021_sub_ice_topo.nc'; 
lat_bathy = ncread(fname, 'lat'); 
lon_bathy = ncread(fname, 'lon'); 

lat_ind = find(lat_bathy >= min(lat)-0.5 & lat_bathy <= max(lat)+0.5); 
lon_ind = find(lon_bathy >= min(lon)-0.5 & lon_bathy <= max(lon)+0.5); 

lat_bathy = lat_bathy(lat_ind); 
lon_bathy = lon_bathy(lon_ind); 

start   = [lon_ind(1) lat_ind(1)]; 
count   = [lon_ind(end) - lon_ind(1) lat_ind(end) - lat_ind(1)] + 1; 
    
bathy_all   = double(ncread(fname, 'elevation', start, count))'; 
        
% now interpolate for the bathymetry transect: 
    % this step flips the sign of the elevation variable. Good if you are
    % plotting depth positive downwards; otherwise it can be easily flipped
    % again afterwards: 
    %   >> bathy = -1 * bathy;
    %   >> yfill = -1 * yfill; 

    bathy = -1*(interp2(lon_bathy, lat_bathy, bathy_all, lon, lat)); 
    
% create xfill and yfill 

    xfill = [x; x(end); x(1); x(1)]; 
    yfill = [bathy; max(bathy); max(bathy); bathy(1)]; 
    
    if sum(isnan(yfill)) > 0 
        yfill = fillmissing(yfill, 'linear'); 
    end
    













end