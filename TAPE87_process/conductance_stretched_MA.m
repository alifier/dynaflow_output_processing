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
%
%   THIS VERSION IS EDITED BY MICHAIL ALIFIERAKIS TO WORK FOR UNSTRUCTURED
%   MESHES

clear all
close all
%RandStream.setDefaultStream(RandStream('mt19937ar','seed',10000));

nofgs = 'n'; % If there are no fgs sheets declare this as 'y'
% Movie?
movieflag = 'yes';
% FILE INPUTS
file = 'S1_dense';
%file = 'Randomtest1';
%file = 'S1';
% Node file
nodefile = strcat(file,'_nodes');
% Elements file
elementfile = strcat(file,'_elements');
% Links file
linksfile = strcat(file,'_springs');
% .dat filename
filename0 = strcat(file,'_trusses');
% TAPE87 extension
filename1 = 's1_dense';
%filename1 = 'Randomtest1';
%filename1 = 's1';
% Side length of mesh (microns)
%l = 5.; 
% Mesh units 
m = 5; % MESH UNIT SEEMS TOO LARGE?
% Number of FGSs
%nb = 1; 
% Number of segments per FGS
lb = 10;
% number of 'steps' i.e. at how many time values have data output from
% Dynaflow? 
%nstep = 22; %increments created by nd plus 2

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


% NODES DEFINITION
node = dlmread(strcat(nodefile,'.dat'));
total_nodes = length(node(:,1));
l =  max(node(:,3)); % length of mesh (microns)

% Elements
elements = dlmread(strcat(elementfile,'.dat'));
elements(:,7) = elements(:,3);
total_elements = length(elements(:,1));

% Load FGS distribution file
if nofgs == 'y'
    beam = [];
    links = [];
    nb = 0;
else
    beam = dlmread(strcat(filename0,'.dat'));
    links = dlmread(strcat(linksfile,'.dat'));
    nb = round(length(beam(:,1))/lb);
end
total_beams = lb*nb;
total_links = (lb+1)*nb;
%total_links = 2*lb*nb; % when lb == 1
%total_links = 2*(lb+1)*nb;

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
%{
al = zeros(1,nstep);
sl = zeros(1,nstep);
conductl = zeros(1,nstep);
conductt = zeros(1,nstep);
conductt3 = zeros(1,nstep);
conductivl = zeros(1,nstep);
conductivt = zeros(1,nstep);
conductivt3 = zeros(1,nstep);
%}

%% READ DISPLACEMENT FILE (TAPE87)

% ktvect = [0:25000:150000]; % w/ tearing
ktvect = []; % Note that this needs to be modified for the case 'with tearing'
ndofs = 2; % Number of dimensions in model
[dm vm rm] = read_data_TAPE87_KS2(['TAPE87.' filename1],'./',ndofs,ktvect,total_nodes);
sizedm=size(dm);
nstep=sizedm(1);

% Use dmm instead of dm if some nodes are not part of the elements
% dmm=zeros(nstep,2*total_nodes);
% dml=length(dm(1,:))/2; %number of nodes that were saved on TAPE87


%lines=[29,46]; %s7, s10, s11, s12 these are the points that the lines start and end
%lines=[29,66]; %s8, s9
%lines=[601,1500]; %Randomtest1

% dmm(:,1:lines(1)-1)=dm(:,1:lines(1)-1);
% dmm(:,lines(2)+1:total_nodes)=dm(:,lines(1):dml);
% 
% dmm(:,total_nodes+1:total_nodes+lines(1)-1)=dm(:,dml+1:dml+lines(1)-1);
% dmm(:,total_nodes+lines(2)+1:end)=dm(:,dml+lines(1):end);
% dm = node displacement vector - first the x displacement of every node is
% listed, then the y displacement, all in one vector
[sigma11,sigma22,sigma33,tau12,eps11,eps22,eps33,gamma12] = read_data_TAPE89(['TAPE89.' filename1],'./',total_elements);
%% LOAD FGS POSITIONS AND CREATE FGS MATRIX

% Define fgs matrix, each row is a sheet with: (start node, end node, start length, end length)
% and node matrix containing all nodes that are part of a sheet 
% stored as (sheet #, all node #'s in individual cols)
if nb>0
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
end
%% Calculate shear stress

trs2=1/3*((sigma11-sigma22).^2+(sigma22-sigma33).^2+(sigma33-sigma11).^2)+2*(tau12.^2);
J2=1/2*trs2;
tau=sqrt(J2);

tre2=1/3*((eps11-eps22).^2+(eps22).^2+(-eps11).^2)+1/2*gamma12.^2;
gamma=sqrt(2*tre2);

mu=tau./gamma; %shear modulus approximation for linear elastic materials
%% LOOP TO TREAT EACH STEP 
while step <= nstep
    
    % Define the positions of the nodes at the step using the displacements
    for i = 1:total_nodes
        nodem(i,1) = i; % same node
        nodem(i,2) = node(i,2) + dm(step,i); % displace in x
        nodem(i,3) = node(i,3) + dm(step,total_nodes+i); % displace in y
%         nodem(i,1) = i; % same node
%         nodem(i,2) = node(i,2) + dmm(step,i); % displace in x
%         nodem(i,3) = node(i,3) + dmm(step,total_nodes+i); % displace in y
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
    %strain(step) = dmm(step,1+total_nodes)*(-2)/l;
    if strain(step) < 1e-3
        strain(step) = 0;
    end
    % width(step) = (l4-l3);

    %% PLOT STRETCHED CONFIGURATION
    % Create new figure and plot the new configuration
        figure(step)
        title(strcat('Strain= ',num2str(strain(step))));
           
        % Plot mesh
        hold on
        for i=1:total_elements
            for j=3:6
                p1=plot([nodem(elements(i,j),2) nodem(elements(i,j+1),2)],[nodem(elements(i,j),3) nodem(elements(i,j+1),3)],'k','Linewidth',0.1);
                p1.Color(4)=0.1;
            end
        end
        
        % Plot teared regions
        if step>1 
            for i=1:total_elements
                %if abs((sigma22(step,i)-sigma22(step-1,i))/(sigma22(step-1,i)))>0.999
                %if sigma22(step,i)/eps22(step,i)<0.001e-3
                if mu(step,i)<3e-8
                    fill([nodem(elements(i,3),2) nodem(elements(i,4),2) nodem(elements(i,5),2) nodem(elements(i,6),2) nodem(elements(i,7),2)],[nodem(elements(i,3),3) nodem(elements(i,4),3) nodem(elements(i,5),3) nodem(elements(i,6),3) nodem(elements(i,7),3)],'r')
                end
            end
        end
        
        % Plot limits of mesh
        x0 = [0 0 l l 0];
        y0 = [0 l l 0 0];
        
        plot(x0,y0,'m','Linewidth',1.5)
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
                plot(x1m,y1m,'b','LineWidth',2.)
                
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
            plot(x1l,y1l,'g','LineWidth',2.) %springs
                
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
      movie2gif(M,[filename1,'.gif'])
end


%% Calculate change in length of the sheets at the end of the stretch
if nb>0
    
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
    
%     length0 = 0;
%     length1 = 0;
%     for i = 1:total_beams
%        length0 = length0 + sqrt((node(beam(i,4),2)-node(beam(i,3),2))^2 + (node(beam(i,4),3)-node(beam(i,3),3))^2);
%        length1 = length1 + sqrt((nodem(beam(i,4),2)-nodem(beam(i,3),2))^2 + (nodem(beam(i,4),3)-nodem(beam(i,3),3))^2);
%        fgs(i,3) = length0;
%        fgs(i,4) = length1;
%     end


    change = fgs(:,4)-fgs(:,3);
    figure
    hist(change);
    set(gca,'FontSize',13);
    xlabel('Change in length of FGS');
    ylabel('Frequency');
    hold on
    % text(0,0,strcat('avg l before: ',num2str(mean(fgs(:,3)),'%6.4g'),' \mum,   \sigma = ',num2str(std(fgs(:,3)),'%6.4g'),' \mum'));
    % text(0,0,strcat('avg l after: ',num2str(mean(fgs(:,4)),'%6.4g'),' \mum,   \sigma = ',num2str(std(fgs(:,4)),'%6.4g'),' \mum'));
end