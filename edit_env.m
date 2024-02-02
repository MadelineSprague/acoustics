% Madeline Sprague (2/1/2024) - mjsprague@uri.edu
%
% This function is nested inside hycom_kr for the unique case that we are
% inputting a particular SSP profile that includes S-wave, density, and
% attenuation parameters (i.e. we include the bottom layers). hycom_kr
% still requires the user to specify a 'bottom' which may be given as the
% bedrock layer underneath the sediment layers specified in full_ssp. 
%
% filename = env filename without file extension as char vector
% full_ssp = double array, with format: 
%   [z(:) cp(:) cs(:) rho(:) ap(:) as(:)] 

function edit_env(filename, full_ssp)

% open .env file and save existing format 

    fid = fopen([filename '.env'], 'r');
    fulltext = strings; 
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        fulltext = [fulltext; string(tline)];
    end
    fulltext = fulltext(2:end); 
    fclose(fid); 

% rewrite certain parts of the text 

    % line 5 - edit bottom depth 
    new_depth = string(full_ssp(end,1)); 
    newstr = split(fulltext(5))'; 
    fulltext(5) = join([newstr(1:2) new_depth]); 

    % change bottom depth in bottom parameters line 
    idx = find(contains(fulltext, "'A'") | contains(fulltext, "''A*''"));  % one line past end of existing SSP 
    newstr = split(fulltext(idx+1));                                       % bottom parameters 
    fulltext(idx+1) = join([new_depth newstr(2:end)']);                    % change bottom depth 
    
    % replace ssp with full_ssp
    fulltext = [fulltext(1:5); join(string(full_ssp),2); fulltext(idx:end)]; % replace ssp with full_ssp 


% rewrite file with old text and with full_ssp 

    fid = fopen([filename '.env'], 'w'); 
    for i = 1:length(fulltext)
        fprintf(fid, [char(fulltext(i)) '\r\n']); 
    end 
    
fclose(fid); 

end