% Madeline Sprague (mjsprague@uri.edu) - October 19 2022
% 
% This function generates a sound speed section using hycom_sec and then
% runs Bellhop using this environmental data and user-set parameters. 
% After running the model it will plot the output. Its usage is the exact 
% same as hycom_sec; after generating a section using hycom_sec, it will 
% put the results through the Bellhop model as desired. Bellhop uses 
% various input and output files, so make sure to set your current 
% directory in a location where you want to store those files. 
%
% USAGE: [outputs] = hycom_bh(lon1, lat1, lon2, lat2, time)
% 
%     After entering the above, the function will generate prompts for 
%     user-set parameters (such as the number of rays, depth of source, 
%     etc.). 
% 
%     Once hycom_bh is complete, it will also produce a struct (outputs)
%     containing all user-set parameters and the data produced by
%     hycom_sec. This allows the user to then re-run the scenario by simply
%     inputting:
%
%       hycom_bh(outputs)
%
%     You may thereby change aspects of the
%     plot by editing outputs directly, instead of re-entering all
%     information via user prompts, or editing .env files - this also
%     reduces runtime by skipping data generation with hycom_sec. For
%     example, after plotting a ray propagation, you may then set 
% 
%       outputs.runtype = '''C'''; 
%       outputs.nrays   = 0; 
%       hycom_bh(outputs);
% 
%     In order to plot the same setting for transmission loss. If you wish
%     to plot a new sound speed field, but retain characteristics of the
%     previous run such as source depth, number of rays, and runtype, you
%     may enter:
%
%       hycom_bh(lat1, lon1, lat2, lon2, time, outputs)
%
%     This will interpolate a new sound speed field from hycom_sec but skip
%     steps requiring manual inputs, all of which are included within
%     outputs. This is the best option for reproducing the same propagation
%     characteristics in different data transects (for example, modeling
%     the same sound source pointing in different directions). 
% 
% Inputs: 
%   lat1, lon1, lat2, lon2: 1x1 double
%       lat1 and lon1 specify starting point of section, lat2 and lon2 
%       specify endpoint. 
%   time: char
%       Specify the time associated with the desired section as a character
%       vector, in a format accepted by datenum: i.e. for August 12th 2021,
%       time = '08-12-2021'. 
%       ***** No time interpolation is needed if you specify a time
%             along a 3-hour interval (00:00:00, 03:00:00, 06:00:00...)
%       ***** Only dates after 12-04-2018 are accepted as a result of HYCOM
%             data availability. 
%   outputs: struct 
%       A struct produced by the function which then may be inputted into
%       the function for another run. If this argument is entered, no other
%       arguments are necessary; however it can be combined with new input
%       coordinates and time as described above. 
%
% Outputs: 
%   outputs: struct 
%       Struct containing all user-set parameters and variables for easy
%       replotting of section. Contains a nested struct also called 
%       'outputs' which stores variables from hycom_sec: 
%           lat, lon: Nx1 double 
%               Arrays giving the lat/lon of each grid point.
%           depth: Mx1 double 
%               Array giving the depth of each grid point.
%           temp: MxN double 
%               Array giving the temperature (deg. C) of each grid point.
%           sal: MxN double 
%               Array giving the salinity (psu) of each grid point.
%           c: MxN double 
%               Array giving the sound speed (m/s) of each grid point.
%           cum_dist: Nx1 double 
%               Array giving the cumulative distance (km) from initial to 
%               final coordinate. 
%           bathy: Nx1 double 
%               Array giving the bathymetry depth (m) of each profile.
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

function [outputs] = hycom_bh(varargin)
    
switch nargin 
    
    case 1 % all input information is contained in outputs struct 
        
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

    % measure to reduce the likelihood of bad extrapolation at depth

        if sum(isnan(outputs.temp(40,:))) > 0 
            outputs.temp(40,:) = mean(outputs.temp(40,:), 'omitnan'); 
            outputs.sal(40,:)  = mean(outputs.sal(40,:),  'omitnan'); 
            outputs.c(40,:)    = mean(outputs.c(40,:),    'omitnan'); 
        end

    % if bathymetry sits below the max depth of data, we need to
    % extrapolate below the bathymetry. 
        
        if max(outputs.bathy) > max(outputs.depth) 
            
            dep_extrap1   = [max(outputs.depth)+50:50:max(outputs.bathy)+50]'; 
            dep_extrap    = repmat(dep_extrap1, 1, width(outputs.c)); 
            pres_extrap   = gsw_p_from_z(-dep_extrap, outputs.lat); 
            sal_extrap    = outputs.sal(end,:);  sal_extrap  = repmat(sal_extrap,  height(dep_extrap), 1); 
            temp_extrap   = outputs.temp(end,:); temp_extrap = repmat(temp_extrap, height(dep_extrap), 1);
            ct_extrap     = gsw_CT_from_t(sal_extrap, temp_extrap, pres_extrap); 
            c_extrap      = gsw_sound_speed(sal_extrap, ct_extrap, pres_extrap); 

            outputs.depth = [outputs.depth; dep_extrap1]; 
            outputs.temp  = [outputs.temp;  temp_extrap]; 
            outputs.sal   = [outputs.sal;   sal_extrap]; 
            outputs.c     = [outputs.c;     c_extrap]; 
            
        end

    % remove all NaNs from sound speed field 

        j = 0; 
        while sum(sum(isnan(outputs.c))) > 0
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
    
        maxdep = max(outputs.depth); 
        ssp    = outputs.c; 
        z      = outputs.depth; 
        bathy  = outputs.bathy; 
        x      = outputs.cum_dist; 
        max_x  = floor(max(outputs.cum_dist)); % seems like max_x needs to be slightly less than where c ends
           if max_x == max(outputs.cum_dist); max_x = max_x - 1; end 
    
    % user-set parameters: 
    
    switch nargin
        
        case 1 % if user has provided inputs 
            
            filename   =  params.filename; 
            env_title  =  params.env_title; 
            runtype    =  params.runtype; 
            freq       =  params.freq; 
            dep_source =  params.dep_source; 
            cp_bot     =  params.cp_bot; 
            cs_bot     =  params.cs_bot; 
            rho_bot    =  params.rho_bot; 
            apt        =  params.apt; 
            ast        =  params.ast; 
            nrays      =  params.nrays; 
            fr_min     = -params.fr_min;
            fr_max     = -params.fr_max;
            n_rec_z    =  params.n_rec_z;
            z_rec_min  =  params.z_rec_min;
            z_rec_max  =  params.z_rec_max;
            n_rec_x    =  params.n_rec_x;
            x_rec_min  =  params.x_rec_min;
            x_rec_max  =  params.x_rec_max;
            
        case 5 % collect manual inputs using a series of prompts
    
            disp('______________________________________________')
            disp('Environmental variables have been successfully generated.') 
            disp('You will now need to input parameters for Bellhop input file creation.')
            disp('All inputs should be entered without apostrophes.')
            disp(' ')

            % general model parameters 
            disp('-----------------')
            disp('General model parameters:')
            filename   = input('Filename (do not include file extension):       ', 's'); 
            env_title  = input('Plot title:                                     ', 's'); 
            env_title  = sprintf('''%s''', env_title); 
            disp('RunType: see options below'); 
                disp('    R = ray trace '); 
                disp('    E = Eigenray trace'); 
                disp('    C = coherent transmission loss');
                disp('    I = incoherent transmission loss'); 
                disp('    A = arrivals'); 
            runtype    = input('Desired RunType:                                ', 's'); 
            runtype    = sprintf('''%s''', runtype); 

            % source parameters 
            disp('-----------------')
            disp('Source parameters:')
            freq       = input('Acoustic frequency (Hz):                        '); 
            dep_source = input('Source depth (m):                               '); 
            nrays      = input('Number of rays (set 0 to let Bellhop choose):   ');
            fr_min     = input('Minimum angle (degrees, unit circle convention: '); 
            fr_max     = input('Maximum angle (degrees, unit circle convention: '); 
            fr_min     = -fr_min; 
            fr_max     = -fr_max; % because BH actually uses a declination convention
            
            % bottom parameters
            disp('-----------------')
            disp('Bottom parameters:')
            cp_bot     = input('P-wave speed (m/s):                             '); 
            cs_bot     = input('S-wave speed (m/s):                             '); 
            rho_bot    = input('Density (kg/m^3):                               '); 
               rho_bot = rho_bot/1000; 
            apt        = input('P-wave attenuation:                             '); 
            ast        = input('S-wave attenuation:                             '); 

            % receiver parameters 
            disp('-----------------')
            disp('Receiver parameters:')
            disp('You can choose to use a grid of receivers spanning the whole')
                disp('sound speed field, or set receiver locations manually.')

            option     = input('Enter ''grid'' or ''manual'':                       ', 's'); 

            switch option 

                case 'grid'

                    n_rec_z    = 100; 
                    z_rec_min  = 0; 
                    z_rec_max  = maxdep;
                    n_rec_x    = 1000;
                    x_rec_min  = 0;
                    x_rec_max  = max_x; 

                case 'manual'

                    n_rec_z    = input('Number of receiver depths:                      '); 
                    z_rec_min  = input('Minimum receiver depth (m):                     '); 
                    z_rec_max  = input('Maximum receiver depth (m):                     '); 
                    n_rec_x    = input('Number of receiver ranges:                      '); 
                    x_rec_min  = input('Minimum receiver range (km):                    '); 
                    x_rec_max  = input('Maximum receiver range (km):                    '); 

            end
            
        case 6 % user has selected a new area but is applying old settings
            
            params     =  varargin{6}; 
            filename   =  params.filename; 
            env_title  =  params.env_title; 
            runtype    =  params.runtype; 
            freq       =  params.freq; 
            dep_source =  params.dep_source; 
            cp_bot     =  params.cp_bot; 
            cs_bot     =  params.cs_bot; 
            rho_bot    =  params.rho_bot; 
            apt        =  params.apt; 
            ast        =  params.ast; 
            nrays      =  params.nrays; 
            fr_min     = -params.fr_min;
            fr_max     = -params.fr_max;
            n_rec_z    =  params.n_rec_z;
            z_rec_min  =  params.z_rec_min;
            z_rec_max  =  params.z_rec_max;
            n_rec_x    =  params.n_rec_x;
            x_rec_min  =  params.x_rec_min;
            x_rec_max  =  params.x_rec_max;

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
    fprintf(fid_env, [char(join(string([maxdep ...
        cp_bot cs_bot rho_bot apt ast]))) ' \r\n']);  % bottom parameters - z, cp, cs, rho, p atten, s atten
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

    formatspec = [repmat('%f ', [1 length(x)]) '\r\n']; 

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

    xfill = [x; x(end); x(1); x(1)]*1000; 
    yfill = [bathy; max(bathy); max(bathy); bathy(1)]; 
        
    % ray plot 

        if strcmp(runtype, '''R''')
            
            figure; hold on

            pcolor(x*1000, z, outputs.c); shading('flat'); colormap(cmocean('speed')); 
            cbar = colorbar; cbar.Label.String = 'Sound Speed (m/s)'; 
            fill(xfill, yfill, [0.8 0.8 0.8])
            plotray([filename '.ray'])
            xlabel('Range (km)')
            ylabel('Depth (m)')
            set(gca, 'ydir', 'rev')
            labelformat([14 16])   
            title(env_title, 'Interpreter', 'none')
            caxis(clims_c); 

            xlim([0 max(xfill)])
            ylim([0 max(yfill)])

            a = gca; 
            xticklabels(string(a.XTick/1000)); 
            
    % arrival plot 

        elseif strcmp(runtype, '''A''')

            plotarr([filename '.arr'], 1, 1, 1) 
            
    % transmission loss plot 

        elseif strcmp(runtype, '''C''') | strcmp(runtype, '''I''')
            
            figure('units', 'normalized', 'position', [0.1 0.1 .8 .8]); hold on 

            plotshd([filename '.shd'])
            fill(xfill, yfill, [.8 .8 .8])
            cbar = colorbar('YDir', 'reverse'); cbar.Label.String = 'Transmission Loss (dB/mkHz)'; 
            colormap(flipud(jet))
            xlabel('Range (km)')
            ylabel('Depth (m)')
            labelformat([14 16])
            title(env_title, 'Interpreter', 'none')
            
            xlim([0 max(xfill)])
            ylim([0 max(yfill)])
            
            a = gca; 
            xticklabels(string(a.XTick/1000)); 

        end
        
% save all data/parameters into outputs struct, so that you can re-run the
% scenario with hycom_bh(outputs)

    outputs = struct('outputs',  outputs,  'lat1',      lat1,      'lon1',       lon1,       ...
                     'lat2',     lat2,     'lon2',      lon2,      'time',       time,       ...
                     'filename', filename, 'env_title', env_title,                           ...
                     'runtype',  runtype,  'freq',      freq,      'dep_source', dep_source, ...
                     'nrays',    nrays,    'fr_min',    -fr_min,   'fr_max',     -fr_max,    ...
                     'n_rec_z',  n_rec_z,  'z_rec_min', z_rec_min, 'z_rec_max',  z_rec_max,  ...
                     'n_rec_x',  n_rec_x,  'x_rec_min', x_rec_min, 'x_rec_max',  x_rec_max,  ...
                     'cp_bot',   cp_bot,   'cs_bot',    cs_bot,    'rho_bot',    rho_bot,    ...
                     'apt',      apt,      'ast',       ast); 
                 
% save outputs into a mat file and return to previous directory 
   
    save([filename '.mat'], 'outputs'); 
    cd(filepath); 

end
