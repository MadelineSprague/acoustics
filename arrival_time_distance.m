% Madeline Sprague (2/1/2024) - mjsprague@uri.edu
%
% This function is somewhat similar to expected_atime, but instead of
% giving the expected arrival time for a given position, it defines the
% arrival-time distance by simulating a line of receivers along the
% transect and finding where the arrival time measurement is along the
% line. This is essentially for the sole case of validating the glider tag
% data. 
%
% USAGE:
% [atime_dist] = arrival_time_distance(lon1, lat1, lon2, lat2, time, atime, sz, rz)
%
% INPUTS:
%   lon1, lat1, lon2, lat2: double
%       coordinates for source (1) and receiver (2)'
%
%   time: datetime or char/string 
%       time pertaining to the transmission for generation of HYCOM data 
%
%   atime: double
%       arrival time measured by receiver in seconds
%
%   sz: double
%       depth of source in meters
%
%   rz: double
%       depth of receiver in meters
%
%   varargin: enter 'plot' in order to see plot of arrival times vs range

function [atime_dist] = arrival_time_distance(lon1, lat1, lon2, lat2, time, atime, sz, rz, varargin)

dir = cd; 
cd('C:\Users\matts\My Drive\Masters\Data\Bellhop_outputs');

    lon2 = lon1 + (lon2-lon1)*1.5; 
    lat2 = lat1 + (lat2-lat1)*1.5; 

    x_rec_min = 1; 
    x_rec_max = m_lldist([lon1 lon2], [lat1 lat2]) * 1.3/1.5;
    n_rec_x   = length(x_rec_min:x_rec_max); 

    outputs = struct('lat1',      lat1,          'lon1',       lon1,      'lat2',      lat2,  'lon2',   lon2,              ...
                     'time',      datestr(time), 'filename',   'atime_dist_temp',                                          ...
                     'env_title', 'atime_dist',  'runtype',    '''A''',                                                    ...
                     'freq',      260,           'dep_source', sz,        'nrays',     0,     'fr_min', -60, 'fr_max', 60, ...
                     'n_rec_z',   1,             'z_rec_min',  rz,        'z_rec_max', rz,                                 ...
                     'n_rec_x',   n_rec_x,       'x_rec_min',  x_rec_min, 'x_rec_max', x_rec_max,                          ...
                     'cp_bot',    1600,          'cs_bot',     600,       'rho_bot',   1800,  'apt',    0.15, 'ast', 0.15); 

% run model and load results 

    outputs =  hycom_bh(lat1, lon1, lat2, lon2, time, outputs);
    [arr, pos] = read_arrivals_asc([outputs.filename '.arr']); 
    range  = pos.r.r; 
    atimes = NaN(size(range)); 

% calculate atime for each receiver position using mean of arrivals
% weighted by amplitude

    for i = 1:length(atimes)

        weight = abs(double(real(arr(i).A))); 
        delay  = double(arr(i).delay); 
        
        atimes(i) = sum(delay .* weight)/sum(weight); 
        %atimes(i) = min(delay); % lets try this out instead

    end

% interpolate for the distance pertaining to the arrival time measurement 

    idx = ~isnan(range) & ~isnan(atimes); 
    atime_dist = interp1(atimes(idx), range(idx)/1000, atime);

% if plots are requested, plot below 

    if strcmp(varargin{:}, 'plot')

        figure; hold on 

            scatter(range/1000, atimes, 'fill')

            p     = polyfit(range, atimes, 1); 
            pvals = polyval(p, range); 
            plot(range/1000, pvals, 'r')
            xlabel('Range (km)')
            ylabel('Arrival Time')

            title(['Average speed: ' char(string(1/(p(1)))) ' m/s'])
            labelformat([14 16])

    end

rmdir(outputs.filename, 's'); 
cd(dir); 

end