% Madeline Sprague (1-30-2024) - mjsprague@uri.edu
%
% This function generally makes use of the hycom_bh framework to generate a
% relationship between distance from source and arrival time. This will
% create and run Bellhop input files using the arrivals RunType, with a
% line of receivers at a given depth. Hence, if the user has an arrival
% time at a known location, we can test how closely this arrival time
% matches the expectation from Bellhop modeling. 

function [outputs] = hycom_bh_arrs(varargin)

switch nargin 
    
    case 1 % all input information is contained in params struct 
        
        params = varargin{:}; 
        
        lat1 = params.lat1; 
        lon1 = params.lon1; 
        lat2 = params.lat2; 
        lon2 = params.lon2; 
        time = params.time; 
        
        outputs = params.outputs; 
        
    case 5 % user inputs lat/lon and time - must run hycom_sec and prompt user for inputs
        
        lon1 = varargin{1}; 
        lat1 = varargin{2};
        lon2 = varargin{3}; 
        lat2 = varargin{4}; 
        time = varargin{5}; 
        
        outputs = hycom_sec(lon1, lat1, lon2, lat2, time); 
        
    case 6 % must generate data for new area, but use old settings and skip prompts
        
        lon1 = varargin{1}; 
        lat1 = varargin{2};
        lon2 = varargin{3}; 
        lat2 = varargin{4}; 
        time = varargin{5}; 
        
        outputs = hycom_sec(lon1, lat1, lon2, lat2, time);    

end
    
% save colorbar limits before doing any extrapolation 

    clims_c = [min(outputs.c, [], 'all') max(outputs.c, [], 'all')];
    
% now prepare data to be written into .env, .bty, and .ssp files. 

    % if bathymetry sits below the max depth of data, we need to
    % extrapolate below the bathymetry. 
        
        if max(outputs.bathy) > max(outputs.depth) 
            
            dep_extrap    = [max(outputs.depth)+50:50:max(outputs.bathy)+50]'; 
            outputs.depth = [outputs.depth; dep_extrap]; 
            dep_extrap    = repmat(dep_extrap, 1, width(outputs.c)); 
            pres_extrap   = gsw_p_from_z(-dep_extrap, outputs.lat); 
            sal_extrap    = outputs.sal(end,:);  sal_extrap  = repmat(sal_extrap,  height(dep_extrap), 1); 
            temp_extrap   = outputs.temp(end,:); temp_extrap = repmat(temp_extrap, height(dep_extrap), 1);
            
            ct_extrap = gsw_CT_from_t(sal_extrap, temp_extrap, pres_extrap); 
            c_extrap  = gsw_sound_speed(sal_extrap, ct_extrap, pres_extrap); 
            outputs.c = [outputs.c; c_extrap]; 
            
        end
        
    % data-dependent parameters: 
    
        maxdep = max(outputs.depth); 
        ssp    = outputs.c; 
        z      = outputs.depth; 
        bathy  = outputs.bathy; 
        x      = outputs.cum_dist; 
        max_x  = floor(max(outputs.cum_dist)); % seems like max_x needs to 
                                               % be slightly less than where c ends
    
    % user-set parameters: 
    
    switch nargin
        
        case 1 % if user has provided inputs 
            
            filename   =  params.filename; 
            env_title  =  params.env_title; 
            runtype    =  params.runtype; 
            freq       =  params.freq; 
            dep_source =  params.dep_source; 
            nrays      =  params.nrays; 
            fr_min     = -params.fr_min;
            fr_max     = -params.fr_max;
            n_rec_z    =  params.n_rec_z;
            z_rec_min  =  params.z_rec_min;
            z_rec_max  =  params.z_rec_max;
            n_rec_x    =  params.n_rec_x;
            x_rec_min  =  params.x_rec_min;
            x_rec_max  =  params.x_rec_max;
            roi        =  params.roi; 
            aoi        =  params.aoi; 
            
        case 5 % needs manual inputs 
    
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
            runtype = '''A'''; 

            % source parameters 
            disp('_____________')
            disp('Source parameters:')
            freq       = input('Acoustic frequency (Hz):                        '); 
            dep_source = input('Source depth (m):                               '); 
            nrays      = 0;
            fr_min     = input('Minimum angle (degrees, unit circle convention: '); 
            fr_max     = input('Maximum angle (degrees, unit circle convention: '); 
            fr_min     = -fr_min; fr_max = -fr_max; % because BH actually uses a declination convention

            % receiver parameters 
            disp('_____________')
            disp(['You may now enter particular arrival times of interest to generate' char(13) ...
                  'associated ranges, or vice versa. Enter these values with only a space between numbers.' char(13) ...
                  'If there are no special values of interest, enter 0.'])
            disp('Receiver parameters:')
            recdep     = input('Enter desired receiver depth (m):               '); 
            aoi        = input('Arrival times (s) of interest:                  '); 
            roi        = input('Ranges (km) of interest:                        ');
            n_rec_z    = 1; 
            z_rec_min  = recdep; 
            z_rec_max  = recdep;
            n_rec_x    = 1000;
            x_rec_min  = 0;
            x_rec_max  = max_x; 
            
        case 6 % user has selected a new area but is applying old settings
            
            params     =  varargin{6}; 
            filename   =  params.filename; 
            env_title  =  params.env_title; 
            runtype    =  params.runtype; 
            freq       =  params.freq; 
            dep_source =  params.dep_source; 
            nrays      =  params.nrays; 
            fr_min     = -params.fr_min;
            fr_max     = -params.fr_max;

            n_rec_z    =  params.n_rec_z;
            z_rec_min  =  params.z_rec_min;
            z_rec_max  =  params.z_rec_max;
            n_rec_x    =  params.n_rec_x;
            x_rec_min  =  params.x_rec_min;
            x_rec_max  =  params.x_rec_max;
            roi        =  params.roi; 
            aoi        =  params.aoi; 

    end
    
% create folder for output files 

    filepath = cd; 
    mkdir(filename)
    addpath([filepath '\' filename])
    cd([filepath '\' filename])    
    
% write .env file 

    fid_env = fopen([filename '.env'], 'w'); 
    fprintf(fid_env, [env_title '\r\n']);             % title 
    fprintf(fid_env, '%3.1f\r\n', freq);              % ray frequency
    fprintf(fid_env, '%d\r\n', 1);                    % NMEDIA
    fprintf(fid_env, '''QVM''\r\n');                  % interpolation/parameters 
    fprintf(fid_env, '%d %1.1f %4.1f\r\n', ...
            [0 0.0 maxdep]);                          % depth limits
    fprintf(fid_env, '%3.2f %5.4f /\r\n', ...
            [z'; ssp(:,1)']);                         % c profile at x = 0
    fprintf(fid_env, '''A*'' 0.0 \r\n');              % acousto-elastic bottom condition 
    fprintf(fid_env, [char(string(maxdep)) ...
            ' 1600.0 0.0 1.0 /\r\n']);                % bottom depth, sound speed of bottom (arbitrary?)
    fprintf(fid_env, [char(string(1)) '\r\n']);       % number of sound sources 
    fprintf(fid_env, '%3.1f /\r\n', dep_source);      % depth of sound source (approximate) [bathy_x(1)-2]
    fprintf(fid_env, [char(string(n_rec_z)) '\r\n']); % number of receiver depths (z)
    fprintf(fid_env, [char(string(z_rec_min)) ' ' ...
                      char(string(z_rec_max)) ...
                      ' /\r\n']);                     % z range of receivers 
    fprintf(fid_env, [char(string(n_rec_x)) '\r\n']); % number of receiver ranges (x)
    fprintf(fid_env, [char(string(x_rec_min)) ' ' ...
                      char(string(x_rec_max)) ...
                      ' /\r\n']);                     % x range of receivers 
    fprintf(fid_env, [runtype '\r\n']);               % RunType ('R' for ray trace)
    fprintf(fid_env, [char(string(nrays)) '\r\n']);   % Number of rays to plot 
    fprintf(fid_env, [char(string(fr_min)) ' ' ...
                      char(string(fr_max)) ...
                      ' /\r\n']);                     % Angular limits of beam fan   
    fprintf(fid_env, ['0.0 ' char(string(maxdep)) ...
                      ' ' char(string(max_x))]);      % ray step size, z limit, x limit 
    fclose(fid_env); 
    
% write .ssp file 

    formatspec = '%f '; 
    for i = 1:length(x)-1
        formatspec = [formatspec '%f ']; 
    end
    formatspec = [formatspec '\r\n']; 

    fid_ssp = fopen([filename '.ssp'], 'w'); 
    fprintf(fid_ssp, '%d\r\n', length(x)); 
    fprintf(fid_ssp, formatspec, x');
    fprintf(fid_ssp, formatspec, ssp'); 
    fclose(fid_ssp); 
    
% write .bty file 

    fid_bty = fopen([filename '.bty'], 'w'); 
    fprintf(fid_bty, 'L\r\n');
    fprintf(fid_bty, '%d\r\n',          length(bathy)); 
    fprintf(fid_bty, '%3.2f %5.4f\r\n', [x bathy]');
    fclose(fid_bty); 

% now run bellhop 

    bellhop(filename); 
    
% plot results 

    figure('units', 'normalized', 'position', [.2 .2 .6 .6]); hold on 

        [arr, pos] = read_arrivals_asc([filename '.arr']); 
        r = pos.r.r/1000;  % km 
        arr1 = NaN(length(arr),1); % record time of all first arrivals
        for i = 1:length(arr)
            scatter(r(i)*ones(arr(i).Narr,1), arr(i).delay, 'fill'); 
            arr1(i) = min(double(arr(i).delay));
        end

        % plot trendline for first arrivals 

            P = polyfit(r, real(arr1), 1); 
            F = polyval(P, r);
            plot(r,F, '--k')

        ylim([0 max(arr1)])
        grid ON 

        % add info from the trendline

            a = gca; 
            dx = (a.XLim(2) - a.XLim(1))/10; 
            dy = (a.YLim(2) - a.YLim(1))/10; 
            text(a.XLim(1)+dx/3, a.YLim(2) - dy/2, ['Average sound speed = ' char(string(1/(P(1))*1000)) ' m/s'], ...
                 'FontSize', 16);
            text(a.XLim(1)+dx/3, a.YLim(2) - dy, ['Endpoint arrival time = ' char(string(min(arr(end).delay))) ' s'], ...
                 'FontSize', 16);

        % add info pertaining to user-specified range or arrival times of
        % interest 

            if roi ~= 0 % plot arrival times associated with ranges of interest 

                roi_arr = polyval(P,roi);

                for i = 1:length(roi)

                    plot([roi(i) roi(i)], [0 roi_arr(i)], '--r')
                    plot([0 roi(i)], [roi_arr(i) roi_arr(i)], '--r')
                    text(0+dx/2, roi_arr(i)+dy/5, ['Range = ' char(string(roi(i))) ' km, arrival = ' char(string(roi_arr(i))) ' sec'])

                end
            end

            if aoi ~=0 % plot ranges associated with arrival times of interest 

                P2 = polyfit(F, r, 1); 
                aoi_r = polyval(P2, aoi); 

                for i = 1:length(aoi)

                    plot([aoi_r(i) aoi_r(i)], [0 aoi(i)], '--r')
                    plot([0 aoi_r(i)], [aoi(i) aoi(i)], '--r')
                    text(0 + dx/2, aoi(i)+dy/5, ['Arrival = ' char(string(aoi(i))) ' sec, range = ' char(string(aoi_r(i))) ' km'])

                end
            end

        xlabel('Range (km)')
        ylabel('Arrival Time (s)')
        title(env_title(2:end-1))
        labelformat([14 16])

% save all data/parameters into outputs struct, so that you can re-run the
% scenario with hycom_bh(outputs)

    outputs = struct('outputs', outputs, 'polyfit', P, 'rec_x', r, 'atime', F, 'roi', roi, 'aoi', aoi, ... 
                     'lat1', lat1, 'lon1', lon1, ...
                     'lat2', lat2, 'lon2', lon2, 'time', time, ...
                     'filename', filename, 'env_title', env_title, ...
                     'runtype', runtype, 'freq', freq, 'dep_source', dep_source, ...
                     'nrays', nrays, 'fr_min', -fr_min, 'fr_max', -fr_max, ...
                     'n_rec_z', n_rec_z, 'z_rec_min', z_rec_min, 'z_rec_max', z_rec_max, ...
                     'n_rec_x', n_rec_x, 'x_rec_min', x_rec_min, 'x_rec_max', x_rec_max); 
                 
% save outputs into a mat file and return to previous directory 
   
    save([filename '.mat'], 'outputs'); 
    cd(filepath); 

end