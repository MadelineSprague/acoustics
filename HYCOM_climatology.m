% Madeline Sprague (2-2-2024) - mjsprague@uri.edu
%
% This script generates a HYCOM climatology of mean conditions from 2019 to
% the present. Main objective is to create two files averaging conditions
% within a large volume of the western North Atlantic, one for summer
% conditions (June, July, August) and another for winter conditions
% (December, January, February). 

% These files are to be used for Bellhop modeling. 

clear 
clc
close all 

%% section 1 - generate winter and summer data 

% choose summer (1) or winter (2) option 

for opt = [1 2]

% define dates 

    summer_dates = [datenum('06-01-2019'):datenum('08-31-2019') ...
                    datenum('06-01-2020'):datenum('08-31-2020') ...
                    datenum('06-01-2021'):datenum('08-31-2021')] + 0.5; % noon every day

    winter_dates = [datenum('12-01-2019'):datenum('02-28-2020') ...
                    datenum('12-01-2020'):datenum('02-28-2021') ...
                    datenum('12-01-2021'):datenum('02-28-2022')] + 0.5; % noon every day

    if opt == 1 
        dates = summer_dates; 
        filename = 'HYCOM_summer.mat'; % 2019-2021
    elseif opt == 2 
        dates = winter_dates; 
        filename = 'HYCOM_winter.mat'; % 2020-2022
    end

% lat lon bounds 

    lat1 = 27.5; 
    lat2 = 42; 
    lon1 = -81.8; 
    lon2 = -57.2; 

% load in hycom data 

    all_c    = cell(length(dates), 1); % collect all information together to easily average
    all_temp = cell(length(dates), 1);
    all_sal  = cell(length(dates), 1);
    n        = length(dates); 
    startnum = 1; 
    looptime = NaN(n,1); 
    failed   = zeros(size(looptime)); 

    for i = startnum:n

        tic
        j = 0; 

        while j < 5

            try
    
                sec_loop(i) = hycom_sec_3d(lat1, lon1, lat2, lon2, dates(i)); 
                all_c{i}    = sec_loop(i).c; 
                all_temp{i} = sec_loop(i).temp; 
                all_sal{i}  = sec_loop(i).sal; 
                j           = 5; 
    
            catch 

                j = j+1; 
                if j == 5; failed(i) = 1; end
                disp(['Run ' char(string(i)) ' failed.'])

            end

        end

% display progress update

    looptime(i) = toc; 
    
    if i ~= n
        time_elapsed = sum(looptime, 'omitnan'); 
        time_remaining = (n-i)/(i+1 - startnum) * time_elapsed;
        if time_remaining > 60*60 
            disp([char(string(i)) '/' char(string(n)) ' runs complete; there are ' ... 
                  char(string(time_remaining/60^2)) ' hours until completion.']);
        else 
            disp([char(string(i)) '/' char(string(n)) ' runs complete; there are ' ...])
                  char(string(time_remaining/60)) ' minutes until completion.']); 
        end
    else 
        disp([char(string(i)) '/' char(string(n)) ' runs complete'])
    end

    end

% display any failed runs 

    if sum(failed) > 0 

        failed = find(failed == 1); 
        disp('The following runs failed after 6 retries:')
        disp(failed)

    end

% concatenate into a large 4-D matrix and average 

    all_c_cat     = cat(4, all_c{:}); 
    all_temp_cat  = cat(4, all_temp{:}); 
    all_sal_cat   = cat(4, all_sal{:}); 

    all_c_mean    = mean(all_c_cat,    4, 'omitnan'); 
    all_temp_mean = mean(all_temp_cat, 4, 'omitnan'); 
    all_sal_mean  = mean(all_sal_cat,  4, 'omitnan'); 

    disp('Averaging Complete'); 

% fill in missing values - including those located beneath bathymetry
    % fill those values first by 2D averaging and then fillmissing on real
    % grid locations

    for i = 1:40

        c_2d    = sec.c(:,:,i);
        temp_2d = sec.temp(:,:,i); 
        sal_2d  = sec.sal(:,:,i); 

        idx          = sec.depth(:,:,i) > abs(sec.bathy) + 100; % fill in values that are well below bathymetry
        c_2d(idx)    = mean(c_2d, 'all'); 
        temp_2d(idx) = mean(temp_2d, 'all'); 
        sal_2d(idx)  = mean(sal_2d, 'all'); 
        
        j = 0; 
        while j < 5

            c_2d = fillmissing(c_2d, 'linear'); % these are more important so they are filled by interpolation
            temp_2d = fillmissing(temp_2d, 'linear');
            sal_2d = fillmissing(sal_2d, 'linear');
            
            if pcnan(c_2d) > 0; j = j+1; elseif pcnan(c_2d) == 0; j = 5; end

        end

        sec.c(:,:,i)    = c_2d; 
        sec.temp(:,:,i) = temp_2d; 
        sec.sal(:,:,i)  = sal_2d; 

    end

% save into a usable file 

    sec = struct('lat',    sec_loop(1).lat,  'lon',   sec_loop(1).lon,  'depth', sec_loop(1).depth, ...
                  'c',     all_c_mean,       'temp',  all_temp_mean,    'sal',   all_sal_mean,      ...
                  'bathy', sec_loop(1).bathy); 
    save(filename, 'sec'); 

end
