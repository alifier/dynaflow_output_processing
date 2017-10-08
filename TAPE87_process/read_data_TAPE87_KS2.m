function [d,varargout] = read_data_TAPE87_KS2(filename,pathname,ndofs,ktvect,nNodes)
% Author: Anastasia Belotsekovets
% Edited and commented by Kevin Sallah
% Most recent edits made: 04.11.2012

% [pathname filename] is the name of the file 
% ndofs is the number of dimensions in the system (number of coordinates
% for each node)
% nNodes is the number of nodes in the system

% Open TAPE87 file
fid = fopen([pathname filename]);

% If you want you can output the size of the file in bytes
% fseek(fid,0,'eof');
% fsize = ftell(fid);
% fseek(fid,0,'bof');

% We want to tell fscanf later to look for some number of coordinates
% written in floating point form 
rstr = '';
for k = 1:ndofs
    rstr = [rstr ' %e'];
end

% The simpler of two methods? We look 
if isempty(ktvect)
    
    kk = 1;
    
    while 1 % Go through file line-by-line until it reaches end 
        % If it reaches end, lines 37-38 will break this loop
        % If it finds the word displacement, loop at 43 will store data
        
        string = fgetl(fid); % Take a new line from the current position
        
        if string == -1 % If you have reached the end, break the loop
            break;
        end
        
        % After seeing the word displacement in a line, it knows that it
        % can start collecting displacement data values from that point on
        if ~isempty(strfind(string,'DISPLACEMENT'))
            % This adds the first value from the line to ktvect 
            ktvect = [ktvect sscanf(string,'%i',1)];
            % Now - read the file from that point, taking only the floating
            % point numbers (the displacements we want) and put them such
            % that the rows are the nodes and the cols are x and y
            % displacements
            A = fscanf(fid,['%*i ' rstr],[ndofs nNodes])';
            % Assemble the output matrix - within each row, x values go first, then y
            d(kk,:) = A(:);
            % Each row in the output matrix is a different time point
            kk = kk+1; 
        end
        
    end

    if nargout>1
        v = [];
        r = [];
    end    
    
end

% Store nothing for velocity or acceleration (we are not using them
% currently)
if nargout>1
    varargout(1) = {v};
    varargout(2) = {r};
end

if nargout == 4
    varargout(3) = {ktvect};
end

fclose('all'); % close file
