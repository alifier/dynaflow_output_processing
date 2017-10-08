function [TIME A B C D E index1 index2] = read_data_TAPE88_2(filename,pathname,nnodes,nodes,nforces,nstresses,stresses,strains)
% Author: Michail Alifierakis


TIME=[];
A=[];
B=[];
C=[];
D=[];
E=[];
% Open TAPE88 file
fid = fopen([pathname filename]);

for i=1:3
    string = fgets(fid);
end

% Read time
counter=0; % How many lines have to be skiped in case data don't need to be stored
while isempty(strfind(string,'Hyperelastic'))==1
    TIME=[TIME; sscanf(string,'%e',8)];
    string=fgets(fid);
    counter=counter+1;
end
kk=0;
ii=0;
index1=zeros(nnodes,1);
index2=zeros(nstresses,1);
cond=1;
while cond==1
    string = fgets(fid);
    if isempty(strfind(string,'NODE'))==0 & numel(strfind(string,'DISPL.  DOF2'))>0
        for i=1:nnodes
            if numel(strfind(string,[' ' num2str(nodes(i)) '    2    4'])) == 1 % the gaps help to not confuse numbers like: 2 with 12 and 21
                kk=kk+1;
                index1(kk)=nodes(i);
                string=fgets(fid);
                while isempty(strfind(string,'Hyperelastic'))==1                    
                    A=[A; sscanf(string,'%e',8)];
                    string=fgets(fid);
                end
            end
        end
        if numel(strfind(string,[' ' num2str(strains) '    2    4'])) == 1
                string=fgets(fid);
                while isempty(strfind(string,'Hyperelastic'))==1                    
                    E=[E; sscanf(string,'%e',8)];
                    string=fgets(fid);
                end
        end
    elseif isempty(strfind(string,'NODE'))==0 & numel(strfind(string,'DISPL.  DOF1'))>0
        
            if numel(strfind(string,[' ' num2str(stresses(1)) '    1    4'])) == 1 % the gaps help to not confuse numbers like: 2 with 12 and 21
                string=fgets(fid);
                while isempty(strfind(string,'Hyperelastic'))==1                    
                    D=[D; sscanf(string,'%e',8)];
                    string=fgets(fid);
                end
            elseif numel(strfind(string,[' ' num2str(stresses(end)) '    1    4'])) == 1
                string=fgets(fid);
                while isempty(strfind(string,'Hyperelastic'))==1                    
                    D=[D; sscanf(string,'%e',8)];
                    string=fgets(fid);
                end
            end
        
    elseif isempty(strfind(string,'NODE'))==0 & numel(strfind(string,'React.  DOF2'))>0
        for i=1:nstresses
            if numel(strfind(string,[' ' num2str(stresses(i)) '    2   15'])) == 1
                ii=ii+1;
                index2(ii)=stresses(i);
                string=fgets(fid);
                while isempty(strfind(string,'Hyperelastic'))==1                    
                    C=[C; sscanf(string,'%e',8)];
                    string=fgets(fid);
                end
            end
        end
    
    elseif numel(strfind(string,'FORC'))==1 
        string=fgets(fid);
        while isempty(strfind(string,'Hyperelastic'))==1                    
            B=[B; sscanf(string,'%e',8)];
            string=fgets(fid);
        end
    
    % skip the rest
    elseif numel(strfind(string,'NODE'))>0
        while isempty(strfind(string,'Hyperelastic'))==1
            string=fgets(fid);
        end
    end
    if string ==  -1 % If you have reached the end, break the loop
        cond=0;
    end
end

fclose('all'); % close file
