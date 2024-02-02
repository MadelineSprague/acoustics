% [outputs] = hycom_sec_3d(lon1, lat1, lon2, lat2, time) produces a 
% 3D grid of temperature and salinity for a given set of lat/lon
% coordinates and a time since 12-04-2018, as well as a 2D grid of seafloor
% bathymetry. Hydrography is interpolated from the HYCOM 3-hourly hindcast
% product; bathymetry is taken from GEBCO, which requires the local file
% GEBCO_2021_sub_ice_topo.nc. Data is pulled from a HYCOM THREDDS server
% hence requiring internet access (it can sometimes take a while). 
%
% This will have issues if your area has longitude values going from
% negative to positive (crossing prime meridian) or the 180 degree line. 
%
% Inputs: 
%   lon1, lat1, lon2, lat2: 1x1 double
%       Both pairs of lat and lon specify opposite corners of the section.
%   time: char
%       Specify the time associated with the desired section as a character
%       vector, in a format accepted by datenum: i.e. for August 12th 2021,
%       time = '08-12-2021'. 
%       ***** Generally less interpolation is needed if you specify a time
%             along a 3-hour interval (00:00:00, 03:00:00, 06:00:00...)
%       An alternative option exists to pull data from my own averaged 
%       HYCOM 2019-2022 files. Enter either 'summer' or 'winter' to get the
%       corresponding file. This file must be available on your MATLAB
%       path. This version ONLY provides the sound speed field (c). 
% Outputs: one output struct, containing the following variables: 
%   lat, lon: MxNx40 double 
%       Arrays giving the lat/lon of each grid point in the section,
%       following meshgrid format at the resolution characteristic of the
%       HYCOM dataset. 
%   depth: MxNx40 double 
%       Array giving the depth of each grid point within the section. 
%   temp: MxNx40 double 
%       Array giving the temperature (deg. C) of each grid point within the 
%       section. 
%   sal: MxNx40 double 
%       Array giving the salinity (psu) of each grid point within the 
%       section. 
%   c: MxNx40 double 
%       Array giving the sound speed (m/s) of each grid point within the
%       section. 
%   bathy: MxN double 
%       Array giving the bathymetry depth (m) of each profile within the
%       section. 

function [outputs] = hycom_sec_3d(lat1, lon1, lat2, lon2, time)

if contains(time, 'summer') | contains(time, 'winter') % user has selected  climatology

    source = ['HYCOM_' time '.mat'];
    if ~exist(source, 'file'); error('Time input does not map to an available climatology file.'); end
    load(source, 'sec'); 

% generate lat/lon vectors 

        latstep = .08; if lat2 < lat1; latstep = -latstep; end
        lonstep = .08; if lon2 < lon1; lonstep = -lonstep; end

        if abs(lat2 - lat1) > abs(lon2 - lon1) 

            lat = [lat1:latstep:lat2]'; 
            lon = [linspace(lon1, lon2, length(lat))]'; 

        elseif abs(lon2 - lon1) > abs(lat2 - lat1) 

            lon = [lon1:lonstep:lon2]'; 
            lat = [linspace(lat1, lat2, length(lon))]'; 

        elseif abs(lat2 - lat1) == abs(lon2 - lon1)

            lat = [lat1:latstep:lat2]'; 
            lon = [lon1:lonstep:lon2]'; 

        end

% 3D grid coordinates 

    depth = squeeze(sec.depth(1,1,:)); 
    [lon_grid, lat_grid, depth_grid] = meshgrid(lon, lat, depth); 

% interpolate 

    warning('off', 'MATLAB:griddedInterpolant:MeshgridEval3DWarnId'); 
    temp = interpn(sec.lon, sec.lat, sec.depth, sec.temp, lon_grid, lat_grid, depth_grid); 
    sal  = interpn(sec.lon, sec.lat, sec.depth, sec.sal,  lon_grid, lat_grid, depth_grid); 
    c    = interpn(sec.lon, sec.lat, sec.depth, sec.c,    lon_grid, lat_grid, depth_grid); 
    warning('on', 'MATLAB:griddedInterpolant:MeshgridEval3DWarnId'); 

else

% URL from which all variables are read 

    % need to update path to include time up until current day 
    
    d = floor(now - datenum('07-31-2021')); 
    t_lim = 7722 + d*8 - 10; 

    path = ['http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?lat[0:1:4250],lon[0:1:4499],' ...
            'time[0:1:' char(string(t_lim)) '],depth[0:1:39],' ...
            'water_temp[0:1:' char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499],' ...
            'salinity[0:1:' char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499]'];

% first we read lat, lon, time, and depth, which are relatively small

    lat_all  = ncread(path, 'lat'); 
    lon_all  = ncread(path, 'lon'); 
    time_all = datenum('01-01-2000') + ncread(path, 'time')/24; 
    depth    = ncread(path, 'depth'); 
    
% change lon convention to -180 to 180 (instead of 0 to 360) 

    lon_all(lon_all > 180) = lon_all(lon_all > 180) - 360; 
    
% calculate step in lat, lon for choosing optimal rectangle of values 

    diff_lat = diff(lat_all); latstep = abs(mean(diff_lat(diff_lat > 0))); 
    diff_lon = diff(lon_all); lonstep = abs(mean(diff_lon(diff_lon > 0))); 
    
% account for lat2/lon2 being less than lat1/lon1 

    if lat2 < lat1; latstep = -latstep; end
    if lon2 < lon1; lonstep = -lonstep; end 
    
% now constrain the temp/sal to only the rectangular area bounded by the
% two corners defined by user inputs 

    lat_ind  = find(lat_all >= min([lat1 lat2])-latstep & lat_all <= max([lat1 lat2])+latstep); 
    lon_ind  = find(lon_all >= min([lon1 lon2])-lonstep & lon_all <= max([lon1 lon2])+lonstep); 
    lat_rect = lat_all(lat_ind); 
    lon_rect = lon_all(lon_ind); 

% and for interpolation, we only need one or two time steps (the times
% before and after, or on, the user specified time) 

    % check for a perfect match 
    
        time_ind = find(time_all == datenum(time)); % if true, no interpolation needed

        if isempty(time_ind) % interpolation needed 

            time_ind = find(time_all - datenum(time) > 0); 
            time_ind = [time_ind(1)-1; time_ind(1)]; 
            
        end
        
    % now load temp and sal variables - this step may take a while 
  
        start = [lon_ind(1)              lat_ind(1)              1   time_ind(1)]; 
        count = [lon_ind(end)-lon_ind(1) lat_ind(end)-lat_ind(1) Inf time_ind(end)-time_ind(1)] + 1;
        temp_rect = ncread(path, 'water_temp', start, count); 
        sal_rect  = ncread(path, 'salinity',   start, count); 
    
    % interpolate through time if necessary 
    
        if count(4) > 1 

            temp_rect = interpn(lon_rect, lat_rect, depth, time_all(time_ind), temp_rect, ...
                                lon_rect, lat_rect, depth, datenum(time)); 
            sal_rect  = interpn(lon_rect, lat_rect, depth, time_all(time_ind), sal_rect,  ...
                                lon_rect, lat_rect, depth, datenum(time)); 
    
        end
        
% remove any fill values 

    temp_rect(temp_rect < -10) = NaN; 
    sal_rect(sal_rect < -10) = NaN;
    
    disp('Data downloaded')

% get sample coordinates for interpolation 

    [lat, lon, depth_s] = meshgrid(lat_rect, lon_rect, depth); 
    
% query points for interpolation: 

    lon_grid = [lon1:lonstep:lon2]'; 
    lat_grid = [lat1:latstep:lat2]'; 
    [lat_grid, lon_grid, depth_grid] = meshgrid(lat_grid, lon_grid, depth); 
    
% interpolate 
    
    temp = interpn(lon, lat, depth_s, temp_rect, ...
                   lon_grid, lat_grid, depth_grid); 
    sal  = interpn(lon, lat, depth_s, sal_rect, ...
                   lon_grid, lat_grid, depth_grid); 
               
% now, need to calculate sound speed using temp/sal/depth 

    pres = gsw_p_from_z(-depth_grid, lat_grid); 
    ct   = gsw_CT_from_t(sal, temp, pres);
    c    = gsw_sound_speed(sal, ct, pres); 

    disp('Variables interpolated');

end

% calculate bathymetry 

    [bathy, ~, ~] = gebco_bathy_3d(lon_grid(:,:,1), lat_grid(:,:,1)); 

% get rid of NaN values in the sound speed field 

    if sum(isnan(c), 'all') > 0 
        c = fillmissing(c, 'linear'); 
        if pcnan(c) == 0
            disp('NaNs removed from sound speed field.')
            disp('Some of the filled-in values may be below the seafloor!')
        elseif pcnan(c) ~= 0 
            disp('Unable to remove NaNs from sound speed field.'); 
        end
    end

% store variables in the output struct 

    outputs = struct('lat',      lat_grid, 'lon',   lon_grid, 'depth', depth_grid, ...
                     'temp',     temp,     'sal',   sal,      'c',     c, ...
                     'bathy', bathy); 
                 
    disp('Variables stored'); 

 end
