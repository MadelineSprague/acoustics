% Madeline Sprague (2/2) - mjsprague@uri.edu
%
% [bathy, lon_grid, lat_grid] = gebco_bathy_3d(lon, lat, varargin) 
% provides bathymetry over an area spanned by the vectors (or grids) lon
% and lat.
% 
% Inputs: 
%   lat, lon: double array 
%       Vectors of lat and lon values for the desired area.
%       Can also be provided as equal-size grid matrices. 
%   3rd argument: 1x3 double (optional)
%       If the user wishes to simply input lat/lon data and have the
%       function create a grid from the limits of these vectors, then they
%       may enter a 3rd argument in the form of a 1x3 double array,
%       containing: 
%           the number of grid points along the smallest dimension
%           the latitude margin (in deg.)
%           the longitude margin (in deg.) 
% 
% Outputs: 
%   bathy: double array 
%       Vector of bathymetry depths for the transect, positive downwards.
%   lat_grid, lon_grid: double array
%       Arrays of lat/lon values produced by the meshgrid of the
%       user-inputted lat and lon vectors, or simply the same as the user
%       inputs if lat and lon are already specified in meshgrid format.

function [bathy, lon_grid, lat_grid] = gebco_bathy_3d(lon, lat, varargin) 

switch nargin 
    case 2

        if isvector(lat) & isvector(lon) % check lat, lon size 
            
            [lat_grid, lon_grid] = meshgrid(lat, lon); 
            
        elseif ismatrix(lat) & ismatrix(lon) & ... 
               ~isvector(lat) & ~isvector(lon) & ...
               size(lat) == size(lon) % lat and lon are both same-size meshgrid arrays 
         
            lat_grid = lat; 
            lon_grid = lon; 
            
        else
            
            error('Provided lat and lon are not compatible.')
            
        end

    case 3 

        threevals  = double(varargin{1}); 
        npoints    = threevals(1); 
        lat_margin = threevals(2); 
        lon_margin = threevals(3); 
        latlims    = [min(lat, [], 'all') max(lat, [], 'all')] + lat_margin * [-1 1]; 
        lonlims    = [min(lon, [], 'all') max(lon, [], 'all')] + lon_margin * [-1 1]; 
        dist_lat   = m_lldist([0 0], latlims); 
        dist_lon   = m_lldist(lonlims, mean(latlims)*[1 1]); 

        if dist_lat <= dist_lon 
            lat = linspace(latlims(1), latlims(2), npoints)'; 
            nlon = round(npoints*dist_lon/dist_lat); 
            lon  = linspace(lonlims(1), lonlims(2), nlon); 
            [lon_grid, lat_grid] = meshgrid(lon, lat); 
        elseif dist_lon < dist_lat
            lon = linspace(lonlims(1), lonlims(2), npoints); 
            nlat = round(npoints*dist_lat/dist_lon); 
            lat  = linspace(latlims(1), latlims(2), nlat); 
            [lon_grid, lat_grid] = meshgrid(lon, lat); 
        end
end
    
% get all lat lon data for the selected region 
    
    fname     = 'GEBCO_2021_sub_ice_topo.nc'; 
    lat_bathy = ncread(fname, 'lat'); 
    lon_bathy = ncread(fname, 'lon'); 

    lat_ind   = find(lat_bathy >= min(lat, [], 'all')-0.005 & lat_bathy <= max(lat, [], 'all')+0.005); 
    lon_ind   = find(lon_bathy >= min(lon, [], 'all')-0.005 & lon_bathy <= max(lon, [], 'all')+0.005); 

    lat_bathy = lat_bathy(lat_ind); 
    lon_bathy = lon_bathy(lon_ind); 

    start     = [lon_ind(1) lat_ind(1)]; 
    count     = [lon_ind(end) - lon_ind(1) lat_ind(end) - lat_ind(1)] + 1; 

    bathy_all = double(ncread(fname, 'elevation', start, count))'; 

% now interpolate for user-specified lat, lon values 

    bathy = interp2(lon_bathy, lat_bathy, bathy_all, lon_grid, lat_grid); 

 end
