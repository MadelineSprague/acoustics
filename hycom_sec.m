% [outputs] = hycom_sec(lon1, lat1, lon2, lat2, time, varargin) produces a 
% temperature and salinity section for a given set of lat/lon coordinates 
% and a time since 12-04-2018. These sections will automatically be plotted 
% with bathymetry as per the user's varargin input. Hydrography is 
% interpolated from the HYCOM 3-hourly hindcast product; bathymetry is 
% taken from GEBCO, which requires the local file
% GEBCO_2021_sub_ice_topo.nc. Data is pulled from a HYCOM THREDDS server
% hence requiring internet access (it can sometimes take a while). 
%
% This function takes a while to run depending on the size of the transect
% - should expect anywhere from 10 seconds for a transect several
% kilometers long, to a minute for transects that are 1000+ km. 
% 
% Inputs: 
%   lon1, lat1, lon2, lat2: 1x1 double
%       lat1 and lon1 specify starting point of section, lat2 and lon2 
%       specify endpoint. 
%   time: char
%       Specify the time associated with the desired section as a character
%       vector, in a format accepted by datenum: i.e. for August 12th 2021,
%       time = '08-12-2021'. 
%       ***** No time interpolation is needed if you specify a time
%             along a 3-hour interval (00:00:00, 03:00:00, 06:00:00...)
%       An alternative option exists to pull data from my own averaged 
%       HYCOM 2019-2022 files. Enter either 'summer' or 'winter' to get the
%       corresponding file. This file must be available on your MATLAB
%       path. This version ONLY provides the sound speed field (c). 
%   varargin (optional)
%       If you wish to plot any or all of the variables, enter any or all
%       of the letters 'tsc' (for temperature, salinity, and sound speed
%       respectively). i.e. if you only want to plot temp and sal, enter
%       'ts', or 'c' for only sound speed. 
%
% Outputs: one output struct, containing the following variables: 
%   lat, lon: Nx1 double 
%       Arrays giving the lat/lon of each grid point within the section. 
%   depth: Mx1 double 
%       Array giving the depth of each grid point within the section. 
%   temp: MxN double 
%       Array giving the temperature (deg. C) of each grid point within the 
%       section. 
%   sal: MxN double 
%       Array giving the salinity (psu) of each grid point within the 
%       section. 
%   c: MxN double 
%       Array giving the sound speed (m/s) of each grid point within the
%       section. 
%   cum_dist: Nx1 double 
%       Array giving the cumulative distance (km) from initial to final
%       coordinate. 
%   bathy: Nx1 double 
%       Array giving the bathymetry depth (m) of each profile within the
%       section. 

function [outputs] = hycom_sec(lon1, lat1, lon2, lat2, time, varargin)

% two possibilities: time is inputted for climatology data, or for the
% 3-hourly HYCOM product 

if isfloat(time); time = datestr(time); end            % convert datenum value to datestring
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

    % interpolate data 

        depth = squeeze(sec.depth(1,1,:)); 
        [lon_grid, depth_grid] = meshgrid(lon', depth); 
        [lat_grid, depth_grid] = meshgrid(lat', depth); 
        c    = interpn(sec.lon, sec.lat, sec.depth, sec.c,    lon_grid, lat_grid, depth_grid); 
        sal  = interpn(sec.lon, sec.lat, sec.depth, sec.sal,  lon_grid, lat_grid, depth_grid); 
        temp = interpn(sec.lon, sec.lat, sec.depth, sec.temp, lon_grid, lat_grid, depth_grid); 

else % use 3-hourly product from the HYCOM THREDDS server 

    % URL from which all variables are read 
    
        % need to update path to include time up until current day 
        
        d = floor(now - datenum('07-31-2021')); 
        t_lim = 7722 + d*8 - 10; 
    
        source = ['http://tds.hycom.org/thredds/dodsC/GLBy0.08/expt_93.0?lat[0:1:4250],lon[0:1:4499],' ...
                'time[0:1:'       char(string(t_lim)) '],depth[0:1:39],'                               ...
                'water_temp[0:1:' char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499],'                 ...
                'salinity[0:1:'   char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499],'                 ...
                'water_u[0:1:'    char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499],'                 ...
                'water_v[0:1:'    char(string(t_lim)) '][0:1:39][0:1:4250][0:1:4499]']; 
    
    % first we read lat, lon, time, and depth, which are relatively small
    
        lat_all  = ncread(source, 'lat'); 
        lon_all  = ncread(source, 'lon'); 
            lon_all(lon_all > 180) = lon_all(lon_all > 180) - 360; 
        time_all = datenum('01-01-2000') + ncread(source, 'time')/24; 
        depth    = ncread(source, 'depth'); 
        
    % now constrain the temp/sal to only the rectangular area bounded by the
    % two corners defined by user inputs 
    
        lat_ind  = find(lat_all >= min([lat1 lat2])-0.5 & lat_all <= max([lat1 lat2])+0.5); 
        lon_ind  = find(lon_all >= min([lon1 lon2])-0.5 & lon_all <= max([lon1 lon2])+0.5); 
        lat_rect = lat_all(lat_ind); 
        lon_rect = lon_all(lon_ind); 
    
    % and for interpolation, we only need one or two time steps (the times
    % before and after, or on, the user specified time) 
    
        % check for a perfect match; otherwise identify the two times from which we will interpolate 
        
            time_ind = find(time_all == datenum(time)); % if true, no interpolation needed
            if isempty(time_ind) % no perfect match - interpolation needed 
    
                time_ind = find(time_all - datenum(time) > 0); 
                time_ind = [time_ind(1)-1; time_ind(1)]; 
                
            end
            
        % now load temp and sal variables - this step may take a while 
        
            try
    
                start = [lon_ind(1)              lat_ind(1)              1   time_ind(1)]; 
                count = [lon_ind(end)-lon_ind(1) lat_ind(end)-lat_ind(1) Inf time_ind(end)-time_ind(1)] + 1;
                temp_rect = ncread(source, 'water_temp', start, count); 
                sal_rect  = ncread(source, 'salinity',   start, count); 
                u_rect    = ncread(source, 'water_u',    start, count);
                v_rect    = ncread(source, 'water_v',    start, count);
    
            catch e 
    
                if contains(e.message, 'Expected start to be positive.')
                    
                    disp('Error: user may have specified a date before 12-04-2018.');
    
                end
                
            end
        
        % interpolate through time if necessary 
        
            if ndims(temp_rect) == 4 
        
                temp_rect = interpn(lon_rect, lat_rect, depth, time_all(time_ind), temp_rect, ...
                                    lon_rect, lat_rect, depth, datenum(time)); 
                sal_rect  = interpn(lon_rect, lat_rect, depth, time_all(time_ind), sal_rect,  ...
                                    lon_rect, lat_rect, depth, datenum(time)); 
                u_rect    = interpn(lon_rect, lat_rect, depth, time_all(time_ind), u_rect,    ...
                                    lon_rect, lat_rect, depth, datenum(time)); 
                v_rect    = interpn(lon_rect, lat_rect, depth, time_all(time_ind), v_rect,    ...
                                    lon_rect, lat_rect, depth, datenum(time)); 
                                
            end
        
    % remove any fill values 
    
        temp_rect(temp_rect < -10) = NaN; 
        sal_rect(sal_rect < -10)   = NaN;
        u_rect(u_rect < -10)       = NaN; 
        v_rect(v_rect < -10)       = NaN; 
        
        disp('Data downloaded')

    % HYCOM has lat/lon on a uniform spacing grid across the globe, so we will
    % follow that resolution while interpolating 
    
        diff_lat = diff(lat_all); latstep = abs(mean(diff_lat(diff_lat > 0))); 
        diff_lon = diff(lon_all); lonstep = abs(mean(diff_lon(diff_lon > 0))); 
        
        % account for lat2/lon2 being less than lat1/lon1 
        
            if lat2 < lat1; latstep = -latstep; end
            if lon2 < lon1; lonstep = -lonstep; end 
        
        lat = [lat1:latstep:lat2]'; 
        lon = [lon1:lonstep:lon2]'; % these are the lat lon values for the section. 
    
        % need to make sure variables are sized the same - normally it's best
        % to limit the resolution to the predefined grid but this produces
        % issues when either lat or lon is constant. 
        
            if length(lat) > length(lon) 
                
                lon = linspace(lon1, lon2, length(lat))';
                
            elseif length(lon) > length(lat) 
                
                lat = linspace(lat1, lat2, length(lon))';
    
            end
    
    % now interpolate to get variables 
    
        [lon_interp, ~]            = meshgrid(lon, depth); 
        [lat_interp, depth_interp] = meshgrid(lat, depth); 
        
        % first: we need to fill in missing values (to look good on figures,
        % and to avoid interpolation errors) - before extrapolating to get
        % unrealistic values that exist below the bathymetry, save colorbar
        % limits pertaining to *real* data 
        
            clims_temp = [floor(min(temp_rect, [], 'all')) ceil(max(temp_rect, [], 'all'))]; 
            clims_sal  = [floor(min(sal_rect,  [], 'all')) ceil(max(sal_rect,  [], 'all'))]; 
        
            temp_rect = fillmissing(temp_rect, 'linear');   
            sal_rect  = fillmissing(sal_rect,  'linear'); 
        
        % we may want a meridional/zonal section where either lat or lon is
        % constant while the other varies. To avoid interpolation errors, do a
        % 1D interp for this case:  
        
            if length(lat_rect) == 1
    
                temp = interp2(lon_rect', depth, squeeze(temp_rect)', lon_interp, depth_interp); 
                sal  = interp2(lon_rect', depth, squeeze(sal_rect)',  lon_interp, depth_interp); 
                u    = interp2(lon_rect', depth, squeeze(u_rect)',    lon_interp, depth_interp); 
                v    = interp2(lon_rect', depth, squeeze(v_rect)',    lon_interp, depth_interp); 
                
            elseif length(lon_rect) == 1
                
                temp = interp2(lat_rect', depth, squeeze(temp_rect)', lat_interp, depth_interp); 
                sal  = interp2(lat_rect', depth, squeeze(sal_rect)',  lat_interp, depth_interp); 
                u    = interp2(lat_rect', depth, squeeze(u_rect)',    lat_interp, depth_interp); 
                v    = interp2(lat_rect', depth, squeeze(v_rect)',    lat_interp, depth_interp); 
    
            else
    
        % or interpolate along varying lat and lon 
    
                temp = interpn(lon_rect, lat_rect, depth, temp_rect, ... 
                               lon_interp, lat_interp, depth_interp, 'linear', -999); 
                           
                sal  = interpn(lon_rect, lat_rect, depth, sal_rect, ...
                               lon_interp, lat_interp, depth_interp, 'linear', -999); 

                u  = interpn(lon_rect, lat_rect, depth, u_rect,     ...
                               lon_interp, lat_interp, depth_interp, 'linear', -999); 

                v  = interpn(lon_rect, lat_rect, depth, v_rect,     ...
                               lon_interp, lat_interp, depth_interp, 'linear', -999); 
                           
            end
            
    % now, need to calculate sound speed using temp/sal/depth 
    
        pres = gsw_p_from_z(-depth_interp, lat); 
        ct   = gsw_CT_from_t(sal, temp, pres);
        c    = gsw_sound_speed(sal, ct, pres); 
    
        disp('Variables interpolated');

end
    
    % calculate bathymetry 
    
        cum_dist = [0; cumsum(m_lldist(lon, lat))]; 
        [bathy, xfill, yfill] = gebco_bathy(lat, lon, cum_dist); 
            
    % store variables in the output struct 
    
        outputs = struct('lat',      lat,      'lon',   lon,   'depth', depth, ...
                         'c',        c,        'sal',   sal,   'temp',  temp,  ...
                         'cum_dist', cum_dist, 'bathy', bathy); 
            if exist('u', 'var'); outputs.u = u; end
            if exist('v', 'var'); outputs.v = v; end
                     
        disp('Variables stored'); 

% plot the transect(s) if requested by user 

    if exist('varargin', 'var')

        if contains(varargin, 't') % temperature section 

            figure('units', 'normalized', 'position', [0 0.4 1/3 0.5])
            hold on 

                colormap(cmocean('thermal')); 
                pcolor(cum_dist, depth, temp); shading('flat');
                colorbar
                labelformat([14 16])
                xlabel('Cumulative Distance (km)')
                ylabel('Depth (m)')
                set(gca, 'ydir', 'rev')
                title(['HYCOM Temperature (' char(176) 'C): ' time])
                caxis(clims_temp); 
                fill(xfill, yfill, [0.8 0.8 0.8])
                xlim([0 max(xfill)])
                ylim([0 max(yfill)])
                
        end

        if contains(varargin, 's') % salinity section 

            figure('units', 'normalized', 'position', [1/3 0.4 1/3 0.5])
            hold on 

                colormap(cmocean('haline'));
                pcolor(cum_dist, depth, sal); shading('flat');
                colorbar
                xlabel('Cumulative Distance (km)')
                ylabel('Depth (m)')
                set(gca, 'ydir', 'rev')
                title(['HYCOM Salinity (psu): ' time])
                labelformat([14 16])       
                caxis(clims_sal); 
                fill(xfill, yfill, [0.8 0.8 0.8])
                xlim([0 max(xfill)])
                ylim([0 max(yfill)])
                
        end
                
        if contains(varargin, 'c') % sound speed section 
            
            % unlike other variables, we need to calculate colorbar limits
            % after extrapolation. We'll do this by taking max and min of
            % all data that resides above the seafloor
            
                ind_clims = zeros(size(depth_interp)); 
                for j = 1:width(depth_interp)
                    
                    ind_clims(:,j) = depth_interp(:,j) < bathy(j); 
                    
                end
                
                ind_clims = logical(ind_clims); 
                clims_c = [floor(min(c(ind_clims), [], 'all')) ceil(max(c(ind_clims), [], 'all'))]; 
            
            figure('units', 'normalized', 'position', [2/3 0.4 1/3 0.5])
            hold on 

                colormap(cmocean('speed'));
                pcolor(cum_dist, depth, c); shading('flat');
                colorbar
                xlabel('Cumulative Distance (km)')
                ylabel('Depth (m)')
                set(gca, 'ydir', 'rev')
                title(['HYCOM Sound Speed (m/s): ' time])
                labelformat([14 16])       
                caxis(clims_c); 
                fill(xfill, yfill, [0.8 0.8 0.8])
                xlim([0 max(xfill)])
                ylim([0 max(yfill)])
                
        end
        
    end

end