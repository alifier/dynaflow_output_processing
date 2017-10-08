% Author: Michail Alifierakis

clear all
%close all
%clc
format long
%%
% Side length of mesh (microns)
%l = 9.4814814815;
%l = 10;
%l = 7.5851851852;
%l=5;
l=9;
% Mesh units 
%m = 100; % MESH UNIT SEEMS TOO LARGE?
% Number of FGSs
nb = 1;
% Number of segments per FGS
lb = 1;

filename1 = 'Aggr_m9p48p0_weak3_3';
%filename1 = 'TEM2_4p0_weak3';
%filename1 = 'Structured6_weak3';
%filename1 = 'Randomtest1_500_strong_slow_short';
%filename1 = 's1_dense_2springs_soft';
%filename1 = 'TEMtear1b';
%filename1 = 'neat3fa';
%filename1 = 'ts0';
%filename1 = 's1_dense';

% how many times nd is ns (relaxation)
%nd = 1;
% final time in the units of dynaflow
%tf=1000;
%dt = 0.1;
%%

% Do you want to plot FGS slippage in y-direction? (y/n)
fgs='n';
% which nodes (state all nodes of polymer matrix and after that all
% equivalent nodes of FGS)

nodes = [5 6 55 56]; % Study 1
%nodes = [5 6 51 52]; % Study 2
%nodes = [5 6 170 171]; % Study 3
%nodes = [5 6 50 51]; % Study 4
%nodes = [5 6 7 8 53 54 55 56]; % Study 5
%nodes = [5 6 7 8 51 52 53 54]; % Study 6
%nodes = [5 6 7 8 260 261 262 263]; % Study 7
%nodes = [5 6 7 8 319 320 321 322]; % Study 8
%nodes = [5 6 7 8 341 342 343 344]; % Study 9
%nodes = [5 6 7 8 249 250 251 252]; % Study 10
%nodes = [5 6 7 8 191 192 193 194]; % Study 11
%nodes = [5 6 7 8 205 206 207 208]; % Study 12
%nodes = [5 6 7 8 388 389 390 391]; % Study 13
%nodes = [5 6 7 8 307 308 309 310]; % Study 14
%nodes = [5 6 7 8 8462 8463 8464 8465]; % Study 15
%nodes = [5 6 7 8 3186 3187 3188 3189]; % Study 16

% Do you want to plot the stresses? (y/n)
stress = 'y';
% which streses to plot (need to put the nodes on top or on bottom, the first and last node need to be the corners)
%stresses = [1,25:33,2]; %TEMt1
%stresses = [1,5:13,2]; %TEMt3
%stresses = [1,7:11,2]; %Studies 1-4
%stresses = [1,9:13,2]; %Studies 5-16
%stresses = [1,205:303,2]; %Randomtest
%stresses = [1,65:83,2]; %Randomtest1_10
%stresses = [1,125:173,2]; %Randomtest1_20
%stresses = [1,245:293,2]; %Randomtest1_40
%stresses = [1,485:533,2]; %Randomtest1_80
%stresses = [1,725:773,2]; %Randomtest1_120 
%stresses = [1,965:1013,2]; %Randomtest1_160
%stresses = [1,1505:1553,2]; %Randomtest1_250
%stresses = [1,3005:3053,2]; %Randomtest1_500
%stresses = [1,413:461,2]; %Structured1
%stresses = [1,311:359,2]; %Structured2
%stresses = [1,167:215,2]; %Structured3
%stresses = [1,617:665,2]; %Structured4
%stresses = [1,35:83,2]; %Structured5
%stresses = [1,443:491,2]; %Structured6
%stresses = [1,419:455,2]; %TEM1_1p0,2p0,4p0
%stresses = [1,647:693,2]; %TEM1_3p0
%stresses = [1,1301:1347,2]; %TEM2_1p0-4p0
%stresses = [1,35:79,2]; %m9p5
%stresses = [1,65:109,2]; %m9p10
%stresses = [1,293:337,2]; %m9p48
%stresses = [1,587:631,2]; %m9p97
%stresses = [1,1175:1219,2]; %m9p195
%stresses = [1,119:163,2]; %m9p19
%stresses = [1,407:451,2]; %m9p67

% Do you want to plot the forces acted by the springs and their displacement? (y/n)
forc = 'n';
% which forces to plot (numbering is according to the number of link element)
%forces = [1 2]; %Studies 1-4
forces = [1 2 3 4]; %Studies 5-16




% How many nodes you need to save 
if fgs == 'y'
    nnodes = numel(nodes);
else
    nnodes = 0;
    nodes = [];
end
% How many lines you want to plot
ndata = 1;



% How many forces you need to save
if forc=='y'
    nforces = numel(forces);
else
    nforces = 0;
    forces = [];
end


if stress == 'y'
    nstresses = numel(stresses);
    strains = stresses(2); %arbitrarily I pick one node on which we will calculate reactions (that means it's on top or bottom)
else
    nstresses = 0;
    stresses = [];
end


%%
color=['b' 'r' 'c' 'm' 'y' 'g' 'k' 'b' 'r' 'c' 'm' 'y' 'g' 'k'];




[TIME,A,B1,B2,Bd1,Bd3,C,D,E,index1,index2,index3] = read_data_TAPE88_4(['TAPE88.' filename1],'./',nnodes,nodes,nforces,forces,nstresses,stresses,strains);

% Matrix A: Stores data for displacements in y-direction for slip 
% calculations - index1 gives the equivalent nodes 
% Matrix B1: Stores data for tangential forces of the springs
% Matrix B2: Stores data for normal forces of the springs
% Matrix Bd1: Stores data for tangential displacement of springs
% Matrix Bd3: Stores data for normal displacement of springs
% Matrix C: Stores data for reactions in y-direction for the stress
% calculations - index2 gives the equivalent nodes
% Matrix D: Stores data for displacements in x-direction for the two edges
% for the stress calculations
% Matrix E: Stores displacement in y for one node at bottom or top for
% strain calculations


% test for reading error
if nnodes~=0
compare1 = eq(index1,nodes');
for i = 1:nnodes
    if compare1(i) == 0
    'reading error at node',nodes(i)
    end
end
end
if nstresses~=0
compare2 = eq(index2,stresses');
for i = 1:nstresses
    if compare2(i) == 0
    'reading error at stress',stresses(i)
    end
end
end

if nforces~=0
compare3 = eq(index3,forces');
for i = 1:nforces
    if compare3(i) == 0
    'reading error at force',forces(i)
    end
end
end

if nnodes ~= 0
    figure(1)
    % which nodes to plot

    node=nodes(1:nnodes/2);
    fgs=nodes((nnodes/2+1):nnodes);
    for j=1:nnodes/2

        for i = 1:nnodes
            if nodes(i) == node(j)
                y = A((i-1)*numel(TIME)+1:i*numel(TIME));
            end
        end

        for i=1:nnodes
            if nodes(i) == fgs(j)
                y = y-A((i-1)*numel(TIME)+1:i*numel(TIME));
            end
        end
        plot(TIME,y,color(j),'LineWidth',2)
        hold on
    end

    ylabel('Distance of Nodes')
    xlabel('Time')
    title('Nodal Slippage in y-direction')


    % edges of FGS, calculate stretching in y
    down(1) = nodes(nnodes/2+1);
    up(1) = down(1) + lb;
    for i=2:nb
        down(i) = up(i-1) + 1;
        up(i) = down(i) + lb;
    end

    ydown=zeros(nb,1);
    yup=zeros(nb,1);
    for j = 1:nb
        for i=1:nnodes
            if down(j) == nodes(i)
                ydown(j)=A(i*numel(TIME));
            end
            if up(j) == nodes(i)
                yup(j)=A(i*numel(TIME));
            end
        end
        stretch(j)=-(ydown(j)-yup(j)); % length of sheet is equal to 1
        % stretch(j)=-(ydown(j)-yup(j))/(l/m*lb);


    %    disp(sprintf('%5d %5d %7g','Stretch % of FGS ',j,stretch(j)*100)) % print stretch percentage of sheet in y

    end

    % create the legend
    M=[];
    for i=1:nnodes/2
        if fix(i/10)==0
            M=[M;['node ',num2str(i)]];
        end
    end

    legend(M)
end


% force figure
if nforces ~= 0
    figure(2)
    for j=1:nforces

        yy = B1((forces(j)-1)*numel(TIME)+1:forces(j)*numel(TIME));

        plot(TIME,yy,color(j),'LineWidth',2)
        hold on
    end
    ylabel('Force')
    xlabel('Time')
    title('Spring Tangential Force')
    legend(M)
    
    figure(3)
    for j=1:nforces

        yy = B2((forces(j)-1)*numel(TIME)+1:forces(j)*numel(TIME));

        plot(TIME,yy,color(j),'LineWidth',2)
        hold on
    end
    ylabel('Force')
    xlabel('Time')
    title('Spring Normal Force')
    legend(M)
end

% spring displacement
if nforces ~= 0
    figure(4)
    for j=1:nforces

        yy = Bd1((forces(j)-1)*numel(TIME)+1:forces(j)*numel(TIME));

        plot(TIME,yy,color(j),'LineWidth',2)
        hold on
    end
    ylabel('DELT-Tangent.')
    xlabel('Time')
    title('Spring Tangential Displacement')
    legend(M)
end

if nforces ~= 0
    figure(5)
    for j=1:nforces

        yy = Bd3((forces(j)-1)*numel(TIME)+1:forces(j)*numel(TIME));

        plot(TIME,yy,color(j),'LineWidth',2)
        hold on
    end
    ylabel('DELT-Normal')
    xlabel('Time')
    title('Spring Normal Displacement')
    legend(M)
end

% stress-strain curve
if nstresses ~= 0
    %figure(6)
    figure(1)
    total_force=0;
    for j=1:nstresses
        total_force=total_force+C((j-1)*numel(TIME)+1:j*numel(TIME));
    end
    Area_true=l-abs(D(1:numel(TIME)))-abs(D(numel(TIME)+1:2*numel(TIME)));
    %true_stress=abs(total_force)*10^6./Area_true;       %trues stress
    %eng_stress=abs(total_force)*10^6/l;
    if any(stresses==1)  % if lower boundary then transform to positive
        eng_stress=-(total_force)*10^6/l;
    else
        eng_stress=(total_force)*10^6/l;
    end
    strain=2*abs(E)/l*100;

    plot(strain,eng_stress)
%     plot(strain,eng_stress,color(i))
    hold on
    set(gca,'fontsize',10)
    %title(filename1)
    xlabel('Strain (%)','fontsize',10)
    ylabel('Engineering Stress (MPa)','fontsize',10)
%    axis([0 100 0 1.6])

    figure(7)
    plot(TIME,eng_stress)
    hold on
    set(gca,'fontsize',10)
    %title(filename1)
    xlabel('Time','fontsize',10)
    ylabel('Engineering Stress (MPa)','fontsize',10)

end


