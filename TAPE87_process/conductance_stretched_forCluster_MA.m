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
%   
% THIS VERSION IS EDITED BY MICHAIL ALIFIERAKIS TO WORK FOR UNSTRUCTURED
%   MESHES AND RUN DIRECTLY ON THE COMPUTER CLUSTER

addpath /home1/malifier/tear/fgs/Random_m9

Directories = ...
    {'m9p5p0','m9p5p2','m9p5p5','m9p5p6',...
    'm9p10p0','m9p10p1','m9p10p2','m9p10p4',...
    'm9p19p0','m9p19p2','m9p19p3','m9p19p8',...
    'm9p48p0','m9p48p2','m9p48p3','m9p48p4',...
    'm9p67p1','m9p67p2','m9p67p4','m9p67p5',...
    'm9p97p2','m9p97p3','m9p97p4','m9p97p5',...
    'm9p195p1','m9p195p3','m9p195p4','m9p195p6'};

torn_area_perc_all = [];
number_of_tears_all = [];

for idirectory=1:numel(Directories)
    
    clearvars -except Directories idirectory torn_area_perc_all number_of_tears_all
    
    directory = Directories(idirectory);

    nofgs = 'n'; % If there are no fgs sheets declare this as 'y'

    % FILE INPUTS
    file = strcat('./',directory,'/Random_',directory);


    % Node file
    nodefile = strcat(file,'_nodes');
    % Elements file
    elementfile = strcat(file,'_elements');
    % Links file
    linksfile = strcat(file,'_springs');
    % .dat filename
    filename0 = strcat(file,'_trusses');
    
    % TAPE87 extension
    filename1 = strcat(file,'_weak3_3_final');

    % Side length of mesh (microns)
    %l = 5.; 
    
    % Mesh units 
    m = 9; % MESH UNIT SEEMS TOO LARGE?
    %m = 5;
    
    % Number of FGSs
    %nb = 1; 
    
    % Number of segments per FGS
    lb = 5;
    %lb = 10;


    r0 = 0.5; 
    r = 2*r0; % the max distance beyond which no conductance is considered
    r1 = 0.04; % the characteristic tunneling distance

    % modulus threshold below which to characterize an element as torn
    torn_threshold = 1.e-7;

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

    step = 1;


    %% READ DISPLACEMENT FILE (TAPE87)

    ktvect = []; % Note that this needs to be modified for the case 'with tearing'
    ndofs = 2; % Number of dimensions in model
    [dm vm rm] = read_data_TAPE87_KS2(['TAPE87.' filename1],'./',ndofs,ktvect,total_nodes);
    sizedm=size(dm);
    nstep=sizedm(1);

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

    %% Calculate area of torn vs non-torn elements

    % torn(time_step,element #), 1 if element is torn, 0 if it's not
    torn = mu<torn_threshold;

    % element_area(time_step,element #)
    % calculate the area of each element for each time step
    element_area = zeros(nstep,total_elements);

    for element=1:total_elements
        % node # for each edge
        node_a = elements(element,3);
        node_b = elements(element,4);
        node_c = elements(element,5);
        node_d = elements(element,6);

        % area of a quad element using determinants (shoelace formula)
        % for triangle the second term reduces to zero because node_c==node_d
        for istep=1:nstep
            % coordinates at step
            xa = node(node_a,2) + dm(istep,node_a);
            ya = node(node_a,3) + dm(istep,total_nodes+node_a);
            xb = node(node_b,2) + dm(istep,node_b);
            yb = node(node_b,3) + dm(istep,total_nodes+node_b);
            xc = node(node_c,2) + dm(istep,node_c);
            yc = node(node_c,3) + dm(istep,total_nodes+node_c);
            xd = node(node_d,2) + dm(istep,node_d);
            yd = node(node_d,3) + dm(istep,total_nodes+node_d);

            % add area of 2 triangles that make up the quad element
            element_area(istep,element) = 1/2*abs(det([xa ya 1;xb yb 1; xc yc 1])) + 1/2*abs(det([xc yc 1;xd yd 1; xa ya 1]));

        end
    end

    total_area = sum(element_area,2);
    torn_area = zeros(nstep,1);
    for istep=1:nstep
        torn_area(istep) = sum(element_area(istep,torn(istep,:)));
    end

    % torn_area_perc(time_step): percentage of torn area relative to total area
    torn_area_perc = torn_area./total_area;
    
    torn_area_perc_all = cat(3,torn_area_perc_all, torn_area_perc);

    %% Calculate number of tears

    % mapping of each node (row index) to its shared elements
    % nodes2elements(total_nodes,10)

    nodes2elements = zeros(total_nodes,10);
    n_nodes2elements = zeros(total_nodes,1);
    for element=1:total_elements

        for i=3:6
            node_i = elements(element,i);

            % if element not already recorded for this node then add it
            if any(nodes2elements(node_i,1:n_nodes2elements(node_i))==element) == 0
                n_nodes2elements(node_i) = n_nodes2elements(node_i) + 1;
                nodes2elements(node_i,n_nodes2elements(node_i)) = element;
            end
        end
    end

    % mapping of each element (row index) to its neighboring elements
    % elements2elements(total_elements,20)

    elements2elements = zeros(total_elements,20);
    n_elements2elements = zeros(total_elements,1);
    for element=1:total_elements
        for i=3:6
            node_i = elements(element,i);
            for j=1:n_nodes2elements(node_i)
                element_j = nodes2elements(node_i,j);
                % if element not itself and not already added in row
                if nodes2elements(node_i,j) ~= element && any(elements2elements(element,1:n_elements2elements(element))==element_j) == 0
                    n_elements2elements(element) = n_elements2elements(element) + 1;
                    elements2elements(element,n_elements2elements(element)) = element_j;
                end
            end
        end
    end

    number_of_tears = zeros(nstep,1);
    % group torn elements in clusters
    for istep=1:nstep
        torn_elements = elements(torn(istep,:),1);
        % cluster at which each element belongs
        cluster = 1:numel(torn_elements);
        checked = [];
        for itorn_element=1:numel(torn_elements)
            element = torn_elements(itorn_element);
            for j=1:numel(checked)
                if any(elements2elements(element,1:n_elements2elements(element)) == checked(j))
                    % re-assign clusters
                    cluster(cluster==cluster(j)) = itorn_element;
                end
            end
            checked = [checked element];
        end
        number_of_tears(istep)=numel(unique(cluster));
    end
    
    number_of_tears_all = cat(3,number_of_tears_all,number_of_tears);

end

save('/home1/malifier/tear/fgs/Random_m9/tear_data.mat')

quit
