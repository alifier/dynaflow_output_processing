function [TIME,A,B1,B2,Bd1,Bd3,C,D,E,index1,index2,index3] = read_data_TAPE88_4(filename,pathname,nnodes,nodes,nforces,forces,nstresses,stresses,strains)
% Author: Michail Alifierakis
% Uses fscanf instead of fgets which makes it faster than
% read_data_TAPE88_2


TIME=[];
A=[];
B1=[];
B2=[];
Bd1=[];
Bd3=[];
C=[];
D=[];
E=[];
% Open TAPE88 file
fid = fopen([pathname filename]);

fscanf(fid,'%s',[5]); % '--' counts as 1 word
Ntsteps=fscanf(fid,'%d',[1])+1;
fscanf(fid,'%s',[1]);

% Read time

TIME=fscanf(fid,'%f',[Ntsteps]);

for i=1:3
    string=fgets(fid);
end
    
kk=0;
ii=0;
jj=0;
index1=zeros(nnodes,1);
index2=zeros(nstresses,1);
index3=zeros(nforces,1);
cond=1;
while cond==1
    if isempty(strfind(string,'NODE'))==0 && numel(strfind(string,'DISPL.  DOF2'))>0
        for i=1:nnodes
            if numel(strfind(string,[' ' num2str(nodes(i)) '    2    4'])) == 1 % the gaps help to not confuse numbers like: 2 with 12 and 21
                kk=kk+1;
                index1(kk)=nodes(i);
                
                                    
                A=[A; fscanf(fid,'%f',[Ntsteps])];
                    
                
            end    
        end
        
        if numel(strfind(string,[' ' num2str(strains) '    2    4'])) == 1
            
            E=fscanf(fid,'%f',[Ntsteps]);
                    
        end
        
    elseif isempty(strfind(string,'NODE'))==0 && numel(strfind(string,'DISPL.  DOF1'))>0
        
            if numel(strfind(string,[' ' num2str(stresses(1)) '    1    4'])) == 1 % the gaps help to not confuse numbers like: 2 with 12 and 21
                
                D=[D; fscanf(fid,'%f',[Ntsteps,1])];
                
            elseif numel(strfind(string,[' ' num2str(stresses(end)) '    1    4'])) == 1
                
                                   
                D=[D; fscanf(fid,'%f',[Ntsteps,1])];
               
            end
        
    elseif isempty(strfind(string,'NODE'))==0 && numel(strfind(string,'React.  DOF2'))>0
        for i=1:nstresses
            if numel(strfind(string,[' ' num2str(stresses(i)) '    2   15'])) == 1
                ii=ii+1;
                index2(ii)=stresses(i);
                
                C=[C; fscanf(fid,'%f',[Ntsteps])];
            end
        end
    
    elseif numel(strfind(string,'FORC'))==1 
        for i=1:nforces
            if numel(strfind(string,['FORC  ELMNT       ' num2str(forces(i))])) == 1
                jj=jj+1;
                index3(jj)=forces(i);
                
                B1 =[B1; fscanf(fid,'%f',[Ntsteps])];
            end
        end
        
    elseif numel(strfind(string,'FORN'))==1 
        for i=1:nforces
            if numel(strfind(string,['FORN  ELMNT       ' num2str(forces(i))])) == 1

                
                B2 =[B2; fscanf(fid,'%f',[Ntsteps])];
            end
        end
        
    elseif numel(strfind(string,'DELT'))==1 
        for i=1:nforces
            if numel(strfind(string,['DELT  ELMNT       ' num2str(forces(i))])) == 1
                Bd1 =[Bd1; fscanf(fid,'%f',[Ntsteps])];
            end
        end
        
    elseif numel(strfind(string,'DELN'))==1 
        for i=1:nforces
            if numel(strfind(string,['DELN  ELMNT       ' num2str(forces(i))])) == 1
                Bd3 =[Bd3; fscanf(fid,'%f',[Ntsteps])];
            end
        end
        
    % skip the rest
    elseif numel(strfind(string,'NODE'))>0
        fscanf(fid,'%f',[Ntsteps]);
    end
    
    string=fgets(fid);
    while isempty(strfind(string,'Hyperelastic'))==1
        string=fgets(fid);
    end
    string=fgets(fid);
    
    if str2num(string) == -1 % If you have reached the end, break the loop
        cond=0;
    end
end

fclose('all'); % close file
