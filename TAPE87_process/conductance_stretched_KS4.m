%    NOT FOR CONDUCTANCE CALCULATION
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Conductance computing for FGS (stretched state)         %
%    written by Anastasia Belotserkovets                     %
%    Version of 1.25.2012                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   THIS VERSION EDITED BY KEVIN SALLAH BEGINNING ON 03.12.2012
%   MOST RECENT EDITS: 05.21.2012
%   THIS VERSION EDITED BY MICHAIL ALIFIERAKIS BEGINNING ON 02.26.2016
%   MOST RECENT EDITS: 

clear all
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',10000));

% FILE INPUTS
% Node file
nodefile = 'S1_nodes';
% Links file
linksfile = 'S1_springs';
% .dat filename
filename0 = 'S1_trusses';
% TAPE87 extension
filename1 = 's0c';
% Side length of mesh (microns)
l = 5.; 
% Mesh units 
m = 10; % MESH UNIT SEEMS TOO LARGE?
% Number of FGSs
nb = 1; 
% Number of segments per FGS
lb = 1;
% number of 'steps' i.e. at how many time values have data output from
% Dynaflow? 
nstep = 2;
% Movie?
movieflag = 'no';

% PARAMETERS
% re = elements number 
% rn = nodes number
% r0 = cutoff
% r = radius
% r1 = tunneling length
% resist0 = maximum resistance
% nb1 = total segments number
% kt = total segment nodes number
% node = nodes [number, x, y]
% nodem = modified nodes [number, x, y]
% elem = elements [number, 4 nodes number]
% nstep = steps number
% conductl = longitudinal conductance
% conductt = transversal conductance

r0 = 0.5; 
r = 2*r0; % the max distance beyond which no conductance is considered
r1 = 0.04; % the characteristic tunneling distance
total_beams = lb*nb;
total_links = 2*(lb+1)*nb;
kt = (lb+1)*nb;
x = zeros(m+1,1);
y = zeros(m+1,1);

% COORDINATES DEFINITION
dmt = l/m;
x(1) = 0;
y(1) = 0;
for i = 1:m
    x(i+1)=x(i)+dmt;
    y(i+1)=y(i)+dmt;
end

% NODES DEFINITION
node = dlmread(strcat(nodefile,'.dat'));
total_nodes = length(node(:,1));

%c = 0;
%yc = 1;
%for i = 1:m+1
%    for j = 1:m+1
%        node(j+c,1) = j+c;
%        node(j+c,2) = x(j);
%        node(j+c,3) = y(yc);
%    end
%    c = c+m+1;
%    yc = yc+1;
%end

step = 1;
al = zeros(1,nstep);
sl = zeros(1,nstep);
conductl = zeros(1,nstep);
conductt = zeros(1,nstep);
conductt3 = zeros(1,nstep);
conductivl = zeros(1,nstep);
conductivt = zeros(1,nstep);
conductivt3 = zeros(1,nstep);


%% READ DISPLACEMENT FILE (TAPE87)

% ktvect = [0:25000:150000]; % w/ tearing
ktvect = []; % Note that this needs to be modified for the case 'with tearing'
ndofs = 2; % Number of dimensions in model
[dm vm rm] = read_data_TAPE87_KS2(['TAPE87.' filename1],'./',ndofs,ktvect,total_nodes);
% dm = node displacement vector - first the x displacement of every node is
% listed, then the y displacement, all in one vector


%% LOAD FGS POSITIONS AND CREATE FGS MATRIX

% Load FGS distribution file
beam = dlmread(strcat(filename0,'.dat'));
links = dlmread(strcat(linksfile,'.dat'));

% Define fgs matrix, each row is a sheet with: (start node, end node, start length, end length)
% and node matrix containing all nodes that are part of a sheet 
% stored as (sheet #, all node #'s in individual cols)
fgs = zeros(nb,4);  
nod = zeros(nb,lb+1);
fgs(1,1) = beam(1,3);
nod(1,1) = beam(1,3);

current_sheet = 1; 
current_node = 1; 

for i = 2:total_beams % check each beam
    if beam(i,3)~=beam(i-1,4) % if the beam is part of the next sheet
       current_sheet = current_sheet+1; 
       fgs(current_sheet-1,2) = beam(i-1,4); % store the end node of previous
       fgs(current_sheet,1) = beam(i,3); % and start node of current sheet
       nod(current_sheet-1,lb+1) = beam(i-1,4);
       nod(current_sheet,1) = beam(i,3);
       current_node = 1; % we are on the first node of a new sheet
    else
       current_node = current_node+1; 
       nod(current_sheet,current_node) = beam(i,3);
    end
end
fgs(current_sheet,2) = beam(total_beams,4); % Store data for very last node
nod(current_sheet,lb+1) = beam(total_beams,4);


%% LOOP TO TREAT EACH STEP 
while step <= nstep
    
    % Define the positions of the nodes at the step using the displacements
    for i = 1:total_nodes
        nodem(i,1) = i; % same node
        nodem(i,2) = node(i,2) + dm(step,i); % displace in x
        nodem(i,3) = node(i,3) + dm(step,total_nodes+i); % displace in y
    end
    
    % I SWITCHED THE METHOD OF CALCULATING STRAIN USED HERE
    % l1 = new bottom
    % l2 = new top
    % l1 = nodem(1,3);
    % l2 = nodem(total_nodes,3);
    % l3 = nodem(1,2);
    % l4 = nodem(total_nodes,2);
    % strain(step) = (l2-l1-l)/l;
    strain(step) = dm(step,1+total_nodes)*(-2)/l;
    if strain(step) < 1e-3
        strain(step) = 0;
    end
    % width(step) = (l4-l3);

    %% PLOT STRETCHED CONFIGURATION
    % Create new figure and plot the new configuration
        figure(step)
        title(strcat('Strain= ',num2str(strain(step))));
           
        % Plot limits of mesh
        x0 = [0 0 l l 0];
        y0 = [0 l l 0 0];
        hold on
        plot(x0,y0,'m')
        hold on
        plot(nodem(:,2),nodem(:,3),'+','MarkerSize',1);
        % text(nodem(1:(m+1)*(m+1),2),nodem(1:(m+1)*(m+1),3),num2str(nodem(1:(m+1)*(m+1),1)),'FontSize',5);
        % text(nodem((m+1)*(m+1):length(nodem(:,1)),2),nodem((m+1)*(m+1):length(nodem(:,1)),3),num2str(nodem((m+1)*(m+1):length(nodem(:,1)),1)),'FontSize',5,'Color','m');
             
        current_beam = 1;
        x1m = zeros(1,lb+1);
        y1m = zeros(1,lb+1);
        
        for i=1:total_beams
            
            x1m(current_beam) = nodem(beam(i,3),2);
            y1m(current_beam) = nodem(beam(i,3),3);
            
            if current_beam == lb;
                x1m(current_beam + 1) = nodem(beam(i,4),2);
                y1m(current_beam + 1) = nodem(beam(i,4),3);
                
                hold on
                plot(x1m,y1m,'b','LineWidth',1.5)
                
                current_beam = 0;
                x1m = zeros(1,lb+1);
                y1m = zeros(1,lb+1);
            end
            
            current_beam = current_beam + 1;
        end
  
        for i=1:total_links
            
            x1l(1) = nodem(links(i,3),2);
            y1l(1) = nodem(links(i,3),3);
            x1l(2) = nodem(links(i,4),2);
            y1l(2) = nodem(links(i,4),3);
                
            hold on
%            plot(x1l,y1l,'g','LineWidth',1.5)
                
        end
        
        axis equal tight
        
 
    % Move on to the next step
    if strcmp(movieflag,'yes') == 1
        M(step) = getframe;
        close all
    end
    step = step+1;
    
end

if strcmp(movieflag,'yes') == 1
      movie(M,1,1)
end


%% Calculate change in length of the sheets at the end of the stretch
ka = 1;
kb = 1;
length0 = 0;
length1 = 0;
for i = 2:total_beams
    if beam(i,3)~=beam(i-1,4)
       ka = ka+1; % move on to the next sheet
       length0 = length0 + sqrt((node(beam(i-1,4),2)-node(beam(i-1,3),2))^2 + (node(beam(i-1,4),3)-node(beam(i-1,3),3))^2);
       length1 = length1 + sqrt((nodem(beam(i-1,4),2)-nodem(beam(i-1,3),2))^2 + (nodem(beam(i-1,4),3)-nodem(beam(i-1,3),3))^2);
       fgs(ka-1,3) = length0;
       fgs(ka-1,4) = length1;
       length0 = 0;
       length1 = 0;
       kb = 1;
    else
       kb = kb+1;
       length0 = length0 + sqrt((node(beam(i,3),2)-node(beam(i-1,3),2))^2 + (node(beam(i,3),3)-node(beam(i-1,3),3))^2);
       length1 = length1 + sqrt((nodem(beam(i,3),2)-nodem(beam(i-1,3),2))^2 + (nodem(beam(i,3),3)-nodem(beam(i-1,3),3))^2);
    end
end
length0 = length0 + sqrt((node(beam(i,4),2)-node(beam(i,3),2))^2 + (node(beam(i,4),3)-node(beam(i,3),3))^2);
length1 = length1 + sqrt((nodem(beam(i,4),2)-nodem(beam(i,3),2))^2 + (nodem(beam(i,4),3)-nodem(beam(i,3),3))^2);
fgs(ka,3) = length0;
fgs(ka,4) = length1;

change = fgs(:,4)-fgs(:,3);
figure
hist(change);
set(gca,'FontSize',13);
xlabel('Change in length of FGS');
ylabel('Frequency');
hold on
% text(0,0,strcat('avg l before: ',num2str(mean(fgs(:,3)),'%6.4g'),' \mum,   \sigma = ',num2str(std(fgs(:,3)),'%6.4g'),' \mum'));
% text(0,0,strcat('avg l after: ',num2str(mean(fgs(:,4)),'%6.4g'),' \mum,   \sigma = ',num2str(std(fgs(:,4)),'%6.4g'),' \mum'));