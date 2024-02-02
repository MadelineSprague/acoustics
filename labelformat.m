% Madeline Sprague (2-2-2024) - mjsprague@uri.edu
%
% Function to auto-format plots with a user-inputted text sizes. This will
% automatically resize titles, axes labels, and legends. 
% 
% Inputs
%   sizes: double array
%       Input the desired fontsizes in the format [labelsize titlesize]. 

function labelformat(sizes)

    h = gca; 

    labelsize = sizes(1); 
    titlesize = sizes(2); 
    
    try

        set(gca, 'FontSize', titlesize); 

        h.XAxis.FontSize = labelsize;
        h.YAxis.FontSize = labelsize;
        
    catch e 
        
        if ~contains(e.message, 'Assigning to 2 elements using a simple assignment')
           
            % this try/catch statement is to deal with an error that pops
            % up when you use labelformat on figures with two y-axes. As
            % far as I can tell the function is still successful, so
            % there's no reason to stop the code. 
            
            return
            
        end
        
    end

end