% Madeline Sprague (1-30-2024) - mjsprague@uri.edu
%
% This function works the same as hycom_bh, loading in data from HYCOM for
% use with the Kraken normal mode acoustic model. See hycom_bh
% documentation for usage instructions. 
%
% Exception: Because KRAKEN is better suited for shallow waveguides than
% HYCOM, I added the capability for the user to specify a custom SSP that
% includes bottom sediment layers. Thus the user can input 
%
%   outputs = hycom_kr(outputs, 'custom', full_ssp)
%
% The outputs struct inputted will first need to be generated on a previous
% hycom_kr run, and full_ssp follows a format described in edit_env.m that
% matches the full six-column SSP block as described in the acoustic
% toolbox documentation. 

function [outputs] = hycom_kr(varargin)
    
switch nargin 
    
    case {1,3} % all input information is contained in outputs struct 
                   % if nargin = 3, 2nd input is 'custom', and the .env SSP
                   % will be replaced with the 3rd input before running
                   % KRAKEN. 
        
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

% KRAKEN is range independent, so we must choose a bottom depth value and
% an SSP representative of the whole area. 

    zbot = round(mean(outputs.bathy));    
        disp(['Standard deviation of bathymetry: ' char(string(std(outputs.bathy)/zbot * 100)) '% of total depth.'])
    ssp  = mean(outputs.c, 2, 'omitnan'); 
        disp(['Average standard deviation of SSP: ' char(string(mean(std(outputs.c, 1, 2), 'omitnan'))) ' m/s.'])
    clow  = 1400; 
    chigh = 1700; 
    

% data-dependent parameters: 

    z      = outputs.depth; 
    x      = outputs.cum_dist; 
    max_x  = floor(max(outputs.cum_dist)); % seems like max_x needs to be slightly less than where c ends
       if max_x == max(outputs.cum_dist); max_x = max_x - 1; end 
    
    % remove NaNs from SSP 

        idx = ~isnan(ssp) & z <= zbot; 
        z   = z(idx); 
        ssp = ssp(idx); 

    % add one seafloor value to SSP 

        ssp_bot = interp1(z, ssp, zbot, 'linear', 'extrap'); 
        if isnan(ssp_bot); ssp_bot = ssp(end); end
        
        z   = [z; zbot]; 
        ssp = [ssp; ssp_bot]; 

        % check for a double value 

            [z, idx, ~] = unique(z); 
            ssp         = ssp(idx); 

% user-set parameters: 
    
    switch nargin
        
        case {1,3} % if user has provided inputs 
            
            filename   = params.filename; 
            env_title  = params.env_title; 
            freq       = params.freq; 
            dep_source = params.dep_source; 
            cp_bot     = params.cp_bot; 
            cs_bot     = params.cs_bot; 
            rho_bot    = params.rho_bot; 
            apt        = params.apt; 
            ast        = params.ast; 
            n_rec_z    = params.n_rec_z;
            z_rec_min  = params.z_rec_min;
            z_rec_max  = params.z_rec_max;
            n_rec_x    = params.n_rec_x;
            x_rec_min  = params.x_rec_min;
            x_rec_max  = params.x_rec_max;
            
        case 5 % collect manual inputs using a series of prompts
    
            disp('______________________________________________')
            disp('Environmental variables have been successfully generated.') 
            disp('You will now need to input parameters for Kraken input file creation.')
            disp('All inputs should be entered without apostrophes.')
            disp(' ')

            % general model parameters 
            disp('-----------------')
            disp('General model parameters:')
            filename   = input('Filename (do not include file extension):       ', 's'); 
            env_title  = input('Plot title:                                     ', 's'); 
                env_title  = sprintf('''%s''', env_title); 

            % source parameters 
            disp('-----------------')
            disp('Source parameters:')
            freq       = input('Acoustic frequency (Hz):                        ');
            dep_source = input('Source depth (m):                               ');
            
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

                    n_rec_z    = 1000; 
                    z_rec_min  = 0; 
                    z_rec_max  = zbot;
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
            freq       =  params.freq; 
            dep_source =  params.dep_source; 

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
    fprintf(fid_env, '''CVMT''\r\n');                 % interpolation/parameters 
    fprintf(fid_env, '%d %1.1f %4.1f\r\n', ...
            [0 0.0 zbot]);                            % depth limits
    fprintf(fid_env, '%3.2f %5.4f /\r\n', ...
            [z'; ssp']);                              % c profile at x = 0
    fprintf(fid_env, '''A'' 0.0 \r\n');               % acousto-elastic bottom condition 
    fprintf(fid_env, [char(join(string([zbot ...
        cp_bot cs_bot rho_bot apt ast]))) ' \r\n']);  % bottom parameters - z, cp, cs, rho, p atten, s atten
    fprintf(fid_env, [char(join(...
        string([clow chigh]))) '\r\n']);              % min and max speeds 
    fprintf(fid_env, [char(string(max_x)) '\r\n']);   % max range (km)
    fprintf(fid_env, [char(string(1)) '\r\n']);       % number of sound sources 
    fprintf(fid_env, '%3.1f /\r\n', dep_source);      % depth of sound source (approximate)
    fprintf(fid_env, [char(string(n_rec_z)) '\r\n']); % number of receiver depths (z)
    fprintf(fid_env, [char(string(z_rec_min)) ' ' ...
                      char(string(z_rec_max)) ...
                      ' /\r\n']);                     % z range of receivers 
    fclose(fid_env); 

% write .flp file 

    fid_flp = fopen([filename '.flp'], 'w'); 
    fprintf(fid_flp, [env_title '\r\n']);                   % title 
    fprintf(fid_flp, '''RA'' \r\n');                        % R = point source, A = adiabatic (C = coupled)
    fprintf(fid_flp, '9999 \r\n');                          % number of modes to include
    fprintf(fid_flp, '1 \r\n');                             % number of profiles 
    fprintf(fid_flp, '0 \r\n');                             % range of profiles
    fprintf(fid_flp, [char(string(n_rec_x)) ' \r\n']);      % number of receivers in x 
    fprintf(fid_flp, ['0 ' char(string(max_x)) ' / \r\n']); % range of receivers
    fprintf(fid_flp, '1 \r\n');                             % number of source depths
    fprintf(fid_flp, [char(string(dep_source)) ' \r\n']);   % source depth (m)
    fprintf(fid_flp, [char(string(n_rec_z)) ' \r\n']);      % number of receiver depths 
    fprintf(fid_flp, ['0 ' char(string(zbot)) '/ \r\n']);   % receiver depths (m)
    fprintf(fid_flp, [char(string(n_rec_z)) ' \r\n']);      % receiver range displacements (?) - must equal NRD
    fprintf(fid_flp, ['0 0 / \r\n']);                       % displacements (zero for vertical array)
    fclose(fid_flp); 

% do we have a custom SSP setup? (I know this isn't very optimized, but it
% was easier to implement this way) 

    if nargin == 3
        if strcmp(varargin{2}, 'custom')
            full_ssp = varargin{3}; 
            edit_env(filename, full_ssp); 
        end
    end

% now run kraken 

    krakenc(filename); 
    field(filename); 

% overlay bathymetry - although the model run assumes a flat bottom 

    xfill = [outputs.cum_dist; outputs.cum_dist(end); outputs.cum_dist(1); outputs.cum_dist(1)]*1000; 
    yfill = [outputs.bathy;    max(outputs.bathy);    max(outputs.bathy);  outputs.bathy(1)]; 
    
% plot results 

    figure('units', 'normalized', 'position', [0.1 0.1 .8 .8]); hold on 
    plotshd([filename '.shd.mat'])
    fill(xfill, yfill, [.75 .75 .75], 'FaceAlpha', 0.3)
    cbar = colorbar('YDir', 'reverse'); cbar.Label.String = 'Transmission Loss (dB/mkHz)'; 
    colormap(flipud(jet))
    xlabel('Range (km)')
    ylabel('Depth (m)')
    labelformat([14 16])
    title(env_title, 'Interpreter', 'none')
    xlim([0 max(x*1000)])
    ylim([0 max(yfill)])
    a = gca; 
    xticklabels(string(a.XTick/1000)); 
        
% save all data/parameters into outputs struct, so that you can re-run the
% scenario with hycom_kr(outputs)

    outputs = struct('outputs',  outputs,  'lat1',       lat1,      'lon1',       lon1,       ...
                     'lat2',     lat2,     'lon2',       lon2,      'time',       time,       ...
                     'filename', filename, 'env_title',  env_title,                           ...
                     'freq',     freq,     'dep_source', dep_source,                          ...
                     'cp_bot',   cp_bot,   'cs_bot',     cs_bot,    'rho_bot',    rho_bot,    ... 
                     'apt',      apt,      'ast',        ast,                                 ...
                     'n_rec_z',  n_rec_z,  'z_rec_min',  z_rec_min, 'z_rec_max',  z_rec_max,  ...
                     'n_rec_x',  n_rec_x,  'x_rec_min',  x_rec_min, 'x_rec_max',  x_rec_max); 
                 
% save outputs into a mat file and return to previous directory 
   
    save([filepath '\' filename '\' filename '.mat'], 'outputs'); 
    cd(filepath); 
    fclose('all'); 

end
