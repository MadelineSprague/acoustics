% This is a 3D counterpart to hycom_bh, which allows for generation of
% HYCOM data in 3D and its use with the Bellhop3D acoustic propagation
% algorithm. It first interpolates HYCOM data for a user-specified lat/lon
% region, or this step may be skipped if the user inputs a parameter struct
% produced by a previous run. It uses the sound speed data and bathymetry
% to produce the needed Bellhop3D input files and run the algorithm,
% producing ray propagations, transmission loss, or arrival plots as
% needed. 
%
% USAGE: [params] = hycom_bh_3d(lat1, lon1, lat2, lon2, time)
% 
%   After entering the above, the function will generate prompts for 
%   user-set parameters (such as the number of rays, depth of source, 
%   etc.). 
% 
%   Once hycom_bh_3d is complete, it will also produce a struct (params)
%   containing all user-set parameters and the data produced by
%   hycom_sec_3d. This allows the user to then re-run the scenario by 
%   simply inputting:
%
%     hycom_bh_3d(params)
%
%   You may thereby change aspects of the
%   plot by editing params directly, instead of re-entering all
%   information via user prompts, or editing .env files - this also
%   reduces runtime by skipping data interpolation with hycom_sec_3d. For
%   example, after plotting a ray propagation, you may then simply set 
% 
%     params.runtype = '''C'''; 
%     params.nalpha  = 0; 
%     params.nbeta   = 0; 
%     hycom_bh(params);
% 
%   In order to plot the same setting for transmission loss. 
% 
% 
% Inputs: 
%   lat, lon1, lat2, lon2: 1x1 double
%       lat1 and lon1 specify starting point of section, lat2 and lon2 
%       specify endpoint. 
%   time: char
%       Specify the time associated with the desired section as a character
%       vector, in a format accepted by datenum: i.e. for August 12th 2021,
%       time = '08-12-2021'. 
%       ***** Generally no interpolation is needed if you specify a time
%             along a 3-hour interval (00:00:00, 03:00:00, 06:00:00...)
%       ***** Only dates after 12-04-2018 are accepted as a result of HYCOM
%             data availability. 
%   params: struct 
%       A struct produced by the function which then may be inputted into
%       the function for another run. If this argument is entered, no other
%       arguments are necessary. 
%
% Outputs: 
%   params: struct 
%       Struct containing all user-set parameters and variables for easy
%       replotting of section. Contains a nested struct called 'outputs'
%       which stores variables from hycom_sec: 
% 
%       lat, lon: MxNx40 double 
%           Arrays giving the lat/lon of each grid point in the section,
%           following meshgrid format at the resolution characteristic of the
%           HYCOM dataset. 
%       depth: MxNx40 double 
%           Array giving the depth of each grid point within the section. 
%       temp: MxNx40 double 
%           Array giving the temperature (deg. C) of each grid point within the 
%           section. 
%       sal: MxNx40 double 
%           Array giving the salinity (psu) of each grid point within the 
%           section. 
%       c: MxNx40 double 
%           Array giving the sound speed (m/s) of each grid point within the
%           section. 
%       bathy: MxN double 
%           Array giving the bathymetry depth (m) of each profile within the
%           section. 
%       
%   Bellhop input files: 
%       [filename].env: the environment file containing user-set
%           parameters.
%       [filename].ssp: the 2D sound speed profile containing user-input
%           sound speed information. 
%       [filename].bty: the bathymetry file containing bathymetry at the
%           range of each profile. 
%
%   Bellhop output files: 
%       [filename].ray, .shd, .arr: the file resulting from running the
%           Bellhop model with the desired RunType. 
%       [filename].prt: the printout file produced by Bellhop to record and
%           display data from the input files. 

function [outputs] = hycom_bh_3d(varargin)

% input arguments (and setting current directory if needed)

    switch nargin 

        case 1 % all input information is contained in params struct 

            params = varargin{:}; 
            fnames = fieldnames(params); 

            for i = 1:length(fnames)

                eval([fnames{i} ' = params.' fnames{i} ';']); 

            end
            
            filepath = cd; 

        case 5 % user inputs lat/lon and time - must run hycom_sec and prompt user for inputs

            lat1 = varargin{1}; 
            lon1 = varargin{2};
            lat2 = varargin{3}; 
            lon2 = varargin{4}; 
            time = varargin{5}; 
            
            disp('______________________________________________')
            disp('The Bellhop algorithm has many file outputs.')
            disp('These files will be stored under the filepath:')
            disp('    Current directory > filename > filename.env, filename.bty, etc.')
            disp('You can paste your desired current directory here or choose your default path,')
            disp('which is set within the function itself.') % default filepath is set within the if loop below 

            filepath = input('Paste your directory here or enter ''default'': ', 's'); 

            if strcmp(filepath, 'default')
                filepath = 'C:\Users\matts\My Drive\Masters\Data\Bellhop_outputs'; % this is your default filepath! 
            end

            if strcmp(filepath, 'cd')
                filepath = cd; 
            end

            cd(filepath)

            disp('Now downloading data...')

            outputs = hycom_sec_3d(lat1, lon1, lat2, lon2, time); 

    end
    
% save colorbar limits before doing any extrapolation 

    clims_c = [min(outputs.c, [], 'all') max(outputs.c, [], 'all')];

% now prepare data to be written into .env, .bty, and .ssp files. 

    % if bathymetry sits below the max depth of data, we need to
    % extrapolate below the bathymetry. Otherwise only extrapolate to 5000
    % m depth 
        
        if max(outputs.bathy, [], 'all') > max(outputs.depth, [], 'all') 
            
            dep_extrap    = [max(outputs.depth, [], 'all')+50:50:max(outputs.bathy, [], 'all')+50]'; 
            sz_lat        = size(outputs.lat);
            dep_extrap    = permute(dep_extrap, [2, 3, 1]); 
            dep_extrap    = repmat(dep_extrap, [sz_lat(1), sz_lat(2), 1]); 
            outputs.depth = cat(3, outputs.depth, dep_extrap); 

            pres_extrap   = gsw_p_from_z(-dep_extrap, outputs.lat(:,:,1:size(dep_extrap, 3))); 
            sal_extrap    = outputs.sal(:,:,end);  sal_extrap  = repmat(sal_extrap,  [1,1,size(dep_extrap, 3)]); 
            temp_extrap   = outputs.temp(:,:,end); temp_extrap = repmat(temp_extrap, [1,1,size(dep_extrap, 3)]); 
            
            ct_extrap = gsw_CT_from_t(sal_extrap, temp_extrap, pres_extrap); 
            c_extrap  = gsw_sound_speed(sal_extrap, ct_extrap, pres_extrap); 
            outputs.c = cat(3, outputs.c, c_extrap); 
            
        end
        
    % remove all NaNs from sound speed field 

        % fill in NaN values that lie below the seafloor 
        for j = 1:width(outputs.depth)
            for i = 1:height(outputs.depth)
                idx_subfloor = outputs.depth(i,j,:) > -outputs.bathy(i,j); 
                outputs.c(i,j,idx_subfloor) = 1550; 
            end
        end

        % now attempt to sweep up missing values 
        j = 0; 
        while sum(isnan(outputs.c), 'all')> 0
            outputs.c = fillmissing(outputs.c, 'linear'); 
            j = j+1; 
            if j == 10 
                disp('Sound speed field NaNs cannot be removed via interpolation;')
                disp('Now attempting extrapolation.')

                    outputs.c = fillmissing(outputs.c, 'nearest', 'EndValues', 'extrap'); 
  
                break
            end
        end

        if sum(isnan(outputs.c), 'all') > 0 
            error('Sound speed NaNs cannot be removed'); 
        end


    % data-dependent parameters: 
    
        maxdep = max(outputs.depth, [], 'all'); 
        ssp    = outputs.c; 
        z      = squeeze(outputs.depth(1,1,:)); 
        bathy  = outputs.bathy; 
        
    % user-set parameters: 
    
        % NOTE: 
            % THETA is the bearing of receivers relative to source 
            % BETA is the bearing of ray fan 
            % ALPHA is the elevation/declination of ray fan 
    
    switch nargin
            
        case 5 % needs manual inputs 

            % need to calculate x and y using both lat and lon - y calculation
            % is based only on lat but x calculation is sensitive to both. So
            % use the latitude value nearest the equator for the lon
            % calculation to produce the longest x line 
            
            % if user inputs params, x/y/z are already specified 

                lat_dummy    = NaN(size(outputs.lon(:,1,1))); 
                lat_dummy(:) = min(abs(outputs.lat), [], 'all'); 
                lon_dummy    = zeros(size(outputs.lat(1,:,1))); 

                x = [0; cumsum(m_lldist(outputs.lon(:,1,1), lat_dummy))]; 
                y = [0; cumsum(m_lldist(lon_dummy, outputs.lat(1,:,1)))]; 

            % for now we can just accept that x will be the same regardless of
            % latitude - i.e. not taking longitude stretching into account
            % based on our latitude 
    
            disp('______________________________________________')
            disp('Data has been successfully produced.') 
            disp('You will now need to input parameters for Bellhop input file creation.')
            disp(' ')

            % general model parameters 
            disp('_____________')
            disp('General model parameters:')
            filename   = input('Filename (do not include file extension):       ', 's'); 
            env_title  = input('Plot title:                                     ', 's'); 
            env_title  = sprintf('''%s''', env_title); 
            disp('RunType: see options below'); disp('    R = ray trace '); 
                                                disp('    C = transmission loss');
                                                disp('    A = arrivals'); 
            runtype    = input('Desired RunType:                                ', 's'); 
            ndims      = input('Enter 2 for Nx2D, 3 for full 3D:                ', 's'); 
            runtype    = [runtype '    ' ndims];
            runtype    = sprintf('''%s''', runtype); 

            % source parameters 
            disp('_____________')
            disp('Source parameters:')
            freq       = input('Acoustic frequency (Hz):                        '); 
            nsx        = input('Number of source coordinates in x (lon):        '); 
            
            if nsx == 1
                sx     = input('Longitude of source (deg):                      '); 
                sx     = interp1(outputs.lon(:,1,1), x, sx);    
            elseif nsx > 1
                sx_min = input('Minimum longitude of source (deg):              '); 
                sx_max = input('Maximum longitude of source (deg):              '); 
                sx_min = interp1(outputs.lon(:,1,1), x, sx_min);   
                sx_max = interp1(outputs.lon(:,1,1), x, sx_max);   
            end
            
            nsy        = input('Number of source coordinates in y (lat):        '); 
            
            if nsy == 1
                sy     = input('Latitude of source (deg):                       '); 
                sy     = interp1(outputs.lat(1,:,1), y, sy); 
            elseif nsy > 1
                sy_min = input('Minimum latitude of source (deg):               '); 
                sy_max = input('Maximum latitude of source (deg):               '); 
                sy_min = interp1(outputs.lat(1,:,1), y, sy_min); 
                sy_max = interp1(outputs.lat(1,:,1), y, sy_max); 
            end
            
            nsz        = input('Number of source coordinates in z (depth, m):   '); 
            
            if nsz == 1
                sz     = input('Depth of source (m):                            '); 
            elseif nsz > 1
                sz_min = input('Minimum depth of source (m):                    '); 
                sz_max = input('Maximum depth of source (m):                    '); 
            end

            disp('You can choose to define beam fan directions manually, or set')
                disp('directions automatically corresponding to an omni-directional source.')
                disp('For each option, you will still choose the number of rays.')
                
            option1    = input('Enter ''omni'' or ''manual'':                       ', 's'); 
            
            nalpha     = input('Number of beams in elevation:                   '); 
            nbeta      = input('Number of beams in bearing direction:           '); 
            
            switch option1
                
                case 'omni' 
                    
                    alpha_min  = -89;
                    alpha_max  = 89;
                    beta_min   = 0; 
                    beta_max   = 360; 
                    
                case 'manual'
            
                    alpha_min  = input('Minimum elevation:                              '); 
                    alpha_max  = input('Maximum elevation:                              '); 
                    alpha_min  = -alpha_min; alpha_max = -alpha_max; % because BH actually uses a declination convention
                    beta_min   = input('Minimum bearing angle (deg):                    '); 
                    beta_max   = input('Maximum bearing angle (deg):                    '); 
            
            end
            
            % redefine x,y to be relative to the source at 0,0 
            
                x = x - sx; 
                y = y - sy; 
                
            % need to specify max x/max y for beam propagation relative to source
            % these are in units of distance (km), not x/y location on the grid. 

                max_x_ray = floor(min(abs([x(1) x(end)])));
                max_y_ray = floor(min(abs([y(1) y(end)])));
            
            % receiver parameters 
            disp('_____________')
            disp('Receiver parameters:')
            disp('You can choose to use a grid of receivers spanning the whole')
                disp('sound speed field, or set receiver locations manually.')
                disp('If you''re not sure, use the grid option.')

            option2     = input('Enter ''grid'' or ''manual'':                       ', 's'); 

            switch option2 

                case 'grid'

                    n_rec_z    = 100; 
                    z_rec_min  = 0; 
                    z_rec_max  = maxdep;
                    n_rec_x    = 100;
                    x_rec_min  = 0;
                    x_rec_max  = max_x_ray; 
                    ntheta     = 37; 
                    theta_min  = 0; 
                    theta_max  = 360; 

                case 'manual'

                    n_rec_z    = input('Number of receiver depths:                      '); 
                    z_rec_min  = input('Minimum receiver depth (m):                     '); 
                    z_rec_max  = input('Maximum receiver depth (m):                     '); 
                    n_rec_x    = input('Number of receiver ranges:                      '); 
                    x_rec_min  = input('Minimum receiver range (km):                    '); 
                    x_rec_max  = input('Maximum receiver range (km):                    '); 
                    ntheta     = input('Number of receiver bearings (deg):              '); 
                    theta_min  = input('Minimum receiver bearing (deg):                 '); 
                    theta_max  = input('Maximum receiver bearing (deg):                 '); 

                    max_x_ray = x_rec_max; 
                    max_y_ray = x_rec_max; 

            end
        
    end
    
% create folder for output files 

    mkdir(filename)
    addpath([filepath '\' filename])
    cd([filepath '\' filename])
    
% write .env file 

    fid_env = fopen([filename '.env'], 'w'); 
    
    fprintf(fid_env, [env_title '\r\n']);             % title 
    fprintf(fid_env, '%3.1f\r\n', freq);              % ray frequency
    fprintf(fid_env, '%d\r\n', 1);                    % NMEDIA
    fprintf(fid_env, '''HVF''\r\n');                  % interpolation/surface condition/TL units 
    fprintf(fid_env, '%d %1.1f %4.1f\r\n', ...
            [0 0.0 maxdep]);                          % depth limits
    fprintf(fid_env, '%3.2f /\r\n', ...
            [z']);                                    % depths (no ssp included)
    fprintf(fid_env, '''A'' 0.0 \r\n');               % acousto-elastic bottom condition, and ~ so it knows to look for bathymetry
    fprintf(fid_env, [char(string(maxdep)) ...
            ' 1600.0 0.0 1.8 0.8 /\r\n']);            % bottom depth, sound speed of bottom (arbitrary?)
    fprintf(fid_env, [char(string(nsx)) '\r\n']);     % number of source coords in x
    
    if nsx == 1 
    fprintf(fid_env, [char(string(0)) ' /\r\n']);     % single x location 
    elseif nsx > 1 
    fprintf(fid_env, [char(string(sx_min)) ' ' ...
                      char(string(sx_max)) ' /\r\n']);% source x range 
    end
    
    fprintf(fid_env, [char(string(nsy)) '\r\n']);     % number of source coords in y 
    
    if nsy == 1 
    fprintf(fid_env, [char(string(0)) ' /\r\n']);     % single y location 
    elseif nsy > 1 
    fprintf(fid_env, [char(string(sy_min)) ' ' ...
                      char(string(sy_max)) ' /\r\n']);% source y range 
    end    
    
    fprintf(fid_env, [char(string(nsz)) '\r\n']);     % number of source coords in z 
    
    if nsz == 1 
    fprintf(fid_env, [char(string(sz)) ' /\r\n']);    % single z location 
    elseif nsz > 1 
    fprintf(fid_env, [char(string(sz_min)) ' ' ...
                      char(string(sz_max)) ...
                      ' /\r\n']);                     % source z range 
    end    

    fprintf(fid_env, [char(string(n_rec_z)) '\r\n']); % number of receiver depths (z)
    fprintf(fid_env, [char(string(z_rec_min)) ' ' ...
                      char(string(z_rec_max)) ...
                      ' /\r\n']);                     % z range of receivers 
    fprintf(fid_env, [char(string(n_rec_x)) '\r\n']); % number of receiver ranges (x)
    fprintf(fid_env, [char(string(x_rec_min)) ' ' ...
                      char(string(x_rec_max)) ...
                      ' /\r\n']);                     % x range of receivers 
                  
    fprintf(fid_env, [char(string(ntheta)) '\r\n']);  % number of receiver bearings 
    fprintf(fid_env, [char(string(theta_min)) ' ' ...
                      char(string(theta_max)) ...
                      ' /\r\n']);                     % range of receiver bearings 
    fprintf(fid_env, [runtype '\r\n']);               % RunType ('R' for ray trace)
    fprintf(fid_env, [char(string(nalpha)) '\r\n']);  % number of beams in elevation/declination 
    fprintf(fid_env, [char(string(alpha_min)) ' ' ...
                      char(string(alpha_max)) ...
                      '/ \r\n']);                     % range of beam elevation angles 
    fprintf(fid_env, [char(string(nbeta)) '\r\n']);   % number of beam bearings 
    fprintf(fid_env, [char(string(beta_min)) ' ' ...
                      char(string(beta_max)) ...
                      '/ \r\n']);                     % range of beam bearings      
    fprintf(fid_env, '%1.1f %f %f %4.1f\r\n', ...
                     [0 max_x_ray max_y_ray maxdep]); % step size, x maximum, y maximum, depth maximum 
    fclose(fid_env); 

% write .ssp file 

    formatspec_x = [repmat('%f ', [1 length(x)]) '\r\n']; 
    formatspec_y = [repmat('%f ', [1 length(y)]) '\r\n']; 
    formatspec_z = [repmat('%f ', [1 length(z)]) '\r\n']; 
    
    fid_ssp = fopen([filename '.ssp'], 'w'); 
    
    fprintf(fid_ssp, '%d\r\n', length(x)); 
    fprintf(fid_ssp, formatspec_x, x'); 
    fprintf(fid_ssp, '%d\r\n', length(y)); 
    fprintf(fid_ssp, formatspec_y, y'); 
    fprintf(fid_ssp, '%d\r\n', length(z)); 
    fprintf(fid_ssp, formatspec_z, z'); 
    
    for i = 1:length(z) 
        
        fprintf(fid_ssp, formatspec_x, ssp(:,:,i)'); 
        
    end
    
    fclose(fid_ssp); 
    
% write bathymetry file 

    bty = -outputs.bathy; 

    fid_bty = fopen([filename '.bty'], 'w'); 
    
    fprintf(fid_bty, '''R''\r\n'); 
    fprintf(fid_bty, '%d\r\n', length(x));
    fprintf(fid_bty, '%f %f / \r\n', [x(1) x(end)]); 
    fprintf(fid_bty, '%d\r\n', length(y));
    fprintf(fid_bty, '%f %f / \r\n', [y(1) y(end)]); 
    
    for i = 1:length(y) 
        
        fprintf(fid_bty, formatspec_x, bty(:,i)'); 
        
    end

    fclose(fid_bty); 
    
% run bellhop3D

    bellhop3d(filename); 
    
% return to previous current directory 

    cd(filepath)
        
% save variables - add sx/sy/sz at the bottom 

    outputs = struct('outputs',   outputs,   'filename',  filename,  'env_title', env_title,                      ...
                     'lat1',      lat1,      'lon1',      lon1,                                                   ...
                     'lat2',      lat2,      'lon2',      lon2,      'time',      time,                           ...
                     'runtype',   runtype,   'freq',      freq,      'ndims',     ndims,                          ...
                     'maxdep',    maxdep,    'nsx',       nsx,       'nsy',       nsy,       'nsz',       nsz,    ...                      ...
                     'nalpha',    nalpha,    'nbeta',     nbeta,                                                  ...
                     'alpha_min', alpha_min, 'alpha_max', alpha_max, 'beta_min',  beta_min, 'beta_max', beta_max, ... 
                     'max_x_ray', max_x_ray, 'max_y_ray', max_y_ray,                                              ... 
                     'n_rec_z',   n_rec_z,   'z_rec_min', z_rec_min, 'z_rec_max', z_rec_max,                      ... 
                     'n_rec_x',   n_rec_x,   'x_rec_min', x_rec_min, 'x_rec_max', x_rec_max,                      ... 
                     'ntheta',    ntheta,    'theta_min', theta_min, 'theta_max', theta_max,                      ...
                     'x',         x,         'y',         y,         'z',         z); 
                 
         if nsx == 1 
             outputs.sx = sx; 
         elseif nsx > 1 
            outputs.sx_min = sx_min; 
            outputs.sx_max = sx_max; 
         end
         
         if nsy == 1 
             outputs.sy = sy; 
         elseif nsy > 1 
            outputs.sy_min = sy_min; 
            outputs.sy_max = sy_max; 
         end
         
         if nsz == 1 
             outputs.sz = sz; 
         elseif nsz > 1 
            outputs.sz_min = sz_min; 
            outputs.sz_max = sz_max; 
         end
            
% now plot the resulting data 

    if contains(runtype, 'R') 
        
        % prepare bathymetry for plotting - add edge values with max depth,
        % to fill in the blank space under the bathymetry when 3D plotting
        
            step = mean(abs(diff(x)))/10; % doesn't really matter to use x or y 

            bathy_x = [min(x)-step; x; max(x)+step]*1000; 
            bathy_y = [min(y)-step; y; max(y)+step]*1000; 
            
            bathy        = -outputs.outputs.bathy'; 
            max_bathy    = max(bathy, [], 'all');                               % edge value 
            bathy        = [NaN(height(bathy), 1) bathy NaN(height(bathy), 1)]; % add L/R limit
            bathy(:,1)   = max_bathy; 
            bathy(:,end) = max_bathy;                                           % add edge value to L/R
            bathy        = [NaN(1, width(bathy)); bathy; NaN(1, width(bathy))]; % add top/bottom limit
            bathy(1,:)   = max_bathy; 
            bathy(end,:) = max_bathy;                                           % add edge value to top/bottom 
            
        figure; hold on 
        colormap(cmocean('topo', 'negative')); 
        labelformat([14 16])
        
        % plot data 
        
            surf(bathy_x, bathy_y, bathy); 
            plotray3d([filename '.ray']) 
        
        % set x limits 
        
            xlim([min(bathy_x) max(bathy_x)])
            ylim([min(bathy_y) max(bathy_y)])
            zlim([0 max(-outputs.outputs.bathy, [], 'all')]) 
            
        % plot formatting 

            grid
            view(30, 60) 
            
            nlabels = 6; 
            
            xtickvals = string(linspace(lon1, lon2, nlabels)); 
            ytickvals = string(linspace(lat1, lat2, nlabels)); 
            
            xtick_labels = cell(1, nlabels); 
            ytick_labels = cell(1, nlabels); 
            for i = 1:nlabels
                
                xtick_labels{i} = sprintf('%3.1f', xtickvals(i)); 
                ytick_labels{i} = sprintf('%3.1f', ytickvals(i)); 
                
            end

            xlabel('Longitude (deg.)')
            ylabel('Latitude (deg.')
            
            xticks(linspace(x(1)*1000, x(end)*1000, nlabels)); 
            yticks(linspace(y(1)*1000, y(end)*1000, nlabels)); 
            
            xticklabels(xtick_labels)
            yticklabels(ytick_labels)
        
    elseif contains(runtype, 'A') 
        
        plotarr3d([filename '.arr'], 1, 1, 1, 1)
        
        
        
    elseif contains(runtype, 'C')
        
        disp('As of now there is no way to plot transmission loss in 3D.')
        disp('This is mainly because the plotshd3d function does not work.')
        disp('Recommend changing runtype and rerunning with output params.')
        
    end

% save the outputs
        
    save([filename '.mat'], 'outputs'); 
    cd(filepath); 

end

