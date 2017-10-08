function [sigma11,sigma22,sigma33,tau12,eps11,eps22,eps33,gamma12] = read_data_TAPE89(filename,pathname,nelem)
% Author: Michail Alifierakis

% [pathname filename] is the name of the file 
% ndofs is the number of dimensions in the system (number of coordinates
% for each node)
% nNodes is the number of nodes in the system

% Open TAPE89 file
fid = fopen([pathname filename]);

% If you want you can output the size of the file in bytes
% fseek(fid,0,'eof');
% fsize = ftell(fid);
% fseek(fid,0,'bof');

% We want to tell fscanf later to look for some number of coordinates
% written in floating point form 
rstr = '';
for i =1:2
    for k = 1:4
        rstr = [rstr ' %e'];
    end
    for k = 1:4
        rstr = [rstr ' %*e'];
    end
    rstr = [rstr ' %*c'];
end

sigma11=zeros(1,nelem);
sigma22=zeros(1,nelem);
sigma33=zeros(1,nelem);
tau12=zeros(1,nelem);
eps11=zeros(1,nelem);
eps22=zeros(1,nelem);
eps33=zeros(1,nelem);
gamma12=zeros(1,nelem);

nstep=1;
string = fgetl(fid);
while 1
    if ~isempty(strfind(string,'STEP'))
        for i=1:nelem
            string = fgetl(fid);
            A = fscanf(fid,rstr,[8 1]);
            sigma11(nstep,i) = A(1);
            sigma22(nstep,i) = A(2);
            sigma33(nstep,i) = A(3);
            tau12(nstep,i) = A(4);
            eps11(nstep,i) = A(5);
            eps22(nstep,i) = A(6);
            eps33(nstep,i) = A(7);
            gamma12(nstep,i) = A(8);
            for j=1:3
                string = fgetl(fid);
            end
        end
        nstep=nstep+1;
    end
    string = fgetl(fid);
    if string == -1
        break;
    end
end

fclose('all'); % close file
