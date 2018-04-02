clc;clear;close all;
 
choice = 4;

switch choice
    case 1
        task1();
    case 2
        task2();
    case 3
        task3a();
        task3b();
    case 4
        task4a();
        task4b();
    case 5
        task5();
    case 6
        task6();
    otherwise
        task1(choice);
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function infoEllipse = task1()
    fprintf('Running Task1 ''\n')
    lambda=500; b=0.02; a1=0.05; a2=0.2; p=0.4; W=[0,1;0,1];
    infoEllipse = rBoolEllipse(lambda,b,a1,a2,p,W);
    fprintf('Task1 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task2()
    clc;clear;close all;
    fprintf('Running Task2 ''\n')
    lambda=200; b=0.02; a1=0.05; a2=0.2; p=0.4; W=[0,1;0,1];
    spacing_1=0.01; spacing_2=0.001;
    rzn_first = rBoolEllipse(lambda,a1,a2,p,b,W);
    rzn_second = rBoolEllipse(lambda,a1,a2,p,b,W);
    
    B1 = digitizeEllSys(rzn_first,W,spacing_1);
    [nx,ny] = size(B1);
    B1 = 1 - B1;
    figure(1);
    imagesc((1:ny),(1:nx),B1);
    colormap(gray);

    B2 = digitizeEllSys(rzn_first,W,spacing_2);
    [nx,ny] = size(B2);
    B2 = 1 - B2;
    figure(2);
    imagesc((1:ny),(1:nx),B2);
    colormap(gray);
    % 
    B3 = digitizeEllSys(rzn_second,W,spacing_1);
    [nx,ny] = size(B3);
    B3 = 1 - B3;
    figure(3);
    imagesc((1:ny),(1:nx),B3);
    colormap(gray);
    % 
    B4 = digitizeEllSys(rzn_second,W,spacing_2);
    [nx,ny] = size(B4);
    B4 = 1 - B4;
    figure(4);
    imagesc((1:ny),(1:nx),B4);
    colormap(gray);
    fprintf('Task2 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 3a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task3a()
    fprintf('Running Task3 ''\n')
    % file to generate realisations and get the minkovski functions
    lambda=50; b=0.04; a1=0.2; a2=0.5; p=0.5; W=[0,1;0,1];
    spacing =0.005; m=20; k=40;
    minkovFnTable = zeros(k+1,4,m+1); %3d Matrix to store each realisation.

    %Generating Realisation and storing to the minkovFnTable 
    for i = 1:m
    rzn = rBoolEllipse(lambda,b,a1,a2,p,W);  % Sample Realisation
    B = digitizeEllSys(rzn,W,spacing);  % digitization of the generated Realisation
    minkovFnTable(:,:,i) =estQMinkowskiFcts(B,k,spacing); %  spacing
    end
    
    %Average of the Minkovski Functional
    for i =2:4
        for j = 1:k+1
            avg =0;
            for k =1:m
                temp = double( minkovFnTable(j,i,k));
                avg = temp + avg;
               % minkovFnTable(j:i:21)= minkovFnTable(j:i:21)+ minkovFnTable(j:i:k);
            end
            minkovFnTable(j,i,m+1)= avg/m;
        end
    end
    % plot the image using the script "Q3_plotMinkFun.m"

%%%%%%%%%%%%%%%%%%%%%%%%
    % file to create the plot for undisturbed minkovski fn table
    m = 20;
    figure(5)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,1),'b');
    grid on; xlabel('r'); ylabel('Aa'); 
    title('First Minkowski Function for s=1');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,2,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'r','LineWidth',1.5); 
    % 
    % % To check the validity of the model
    % plot(anal(:,1),anal(:,2),'g');    
    % plot(MiTab(:,1),MiTab(:,2),'y');%                                          
    %                                                             
    hold off
    %%%%%%% 2 Minkowski Function
    figure(6)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,1),'b');
    grid on; xlabel('r'); ylabel('La'); 
    title('Second Minkowski Function for s=1');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,3,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'r','LineWidth',1.5);    
    hold off
    %%%%%%% 3 Minkowski Function
    figure(7)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,1),'b');
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function for s=1');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,4,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'r','LineWidth',1.5);   
    hold off
    fprintf('Task3a Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 3b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task3b()
    % file to generate realisations and get the minkovski functions
    lamda=50; b=0.04; a1=0.2; a2=0.5; p=0.5; W=[0,2;0,2];
    spacing =0.005; m=20; k=40;
    minkovFnTable = zeros(k+1,4,m+1); %3d Matrix to store each realisation.

    %Generating Realisation and storing to the minkovFnTable 
    for i = 1:m
    rzn = rBoolEllipse(lamda,b,a1,a2,p,W);  % Sample Realisation
    B = digitizeEllSys(rzn,W,spacing);  % digitization of the generated Realisation
    minkovFnTable(:,:,i) =estQMinkowskiFcts(B,k,spacing); %  spacing
    end
    
    %Average of the Minkovski Functional
    for i =2:4
        for j = 1:k+1
            avg =0;
            for k =1:m
                temp = double( minkovFnTable(j,i,k));
                avg = temp + avg;
               % minkovFnTable(j:i:21)= minkovFnTable(j:i:21)+ minkovFnTable(j:i:k);
            end
            minkovFnTable(j,i,m+1)= avg/m;
        end
    end
    % plot the image using the script "Q3_plotMinkFun.m"

%%%%%%%%%%%%%%%%%%%%%%%%
    % file to create the plot for undisturbed minkovski fn table
    m = 20;
    figure(8)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,1),'b');
    grid on; xlabel('r'); ylabel('Aa'); 
    title('First Minkowski Function for s=2');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,2,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'r','LineWidth',1.5); 
    % 
    % % To check the validity of the model
    % plot(anal(:,1),anal(:,2),'g');    
    % plot(MiTab(:,1),MiTab(:,2),'y');%                                          
    %                                                             
    hold off
    %%%%%%% 2 Minkowski Function
    figure(9)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,1),'b');
    grid on; xlabel('r'); ylabel('La'); 
    title('Second Minkowski Function for s=2');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,3,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'r','LineWidth',1.5);    
    hold off
    %%%%%%% 3 Minkowski Function
    figure(10)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,1),'b');
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function for s=2');
    hold on
    for i= 2:m
        plot(minkovFnTable(:,1,1),minkovFnTable(:,4,i),'b');
    end
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'r','LineWidth',1.5);   
    hold off
    fprintf('Task3b Completed ''\n')
    fprintf('Task3 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 4a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task4a()
    fprintf('Running Task4 ''\n')
    lambda=50; b=0.04; a1=0.2; a2=0.5; p=0.5; W=[0,1;0,1];
    m=20; k=40;
    spacing=0.005; q=0.001;  % Disturbance coefficient
    minkovFnTable = zeros(k+1,4,m+1);
    distMinkovFnTable = zeros(k+1,4,m+1);

    for i = 1:m    
        rzn = rBoolEllipse(lambda,b,a1,a2,p,W);  % Sample Realisation
        rbool =digitizeEllSys(rzn,W,spacing);  % digitization of the generated Realisation
        minkovFnTable(:,:,i) =estQMinkowskiFcts(rbool,k,spacing); 
        rbool1 = disturbImage(rbool,q); % Disturbing the image
        distMinkovFnTable(:,:,i) =estQMinkowskiFcts(rbool1,k,spacing);
    end

    %minkovFnTable(:,:)
    %Average of the Minkovski Functional
    for i =2:4
        for j = 1:k+1
            avg =0; avg1=0;
            for k =1:m
                temp = double( minkovFnTable(j,i,k));
                temp1 = double( distMinkovFnTable(j,i,k));
                avg = temp + avg;
                avg1 = temp1 + avg1;
               % minkovFnTable(j:i:21)= minkovFnTable(j:i:21)+ minkovFnTable(j:i:k);

            end
            minkovFnTable(j,i,m+1)= avg/m;
            distMinkovFnTable(j,i,m+1)= avg1/m;
        end
    end
    %%%%%%%%%%%%%%%%
    % file to create the plot for undisturbed minkovski fn table
    m = 20;

    figure(11)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,1),'r');
    grid on; xlabel('r'); ylabel('Aa'); 
    title('First Minkowski Function (q=0.001)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,1),'b');
    for i= 2:m

        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,2,i),'b');

    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,m+1),'r','LineWidth',1.5);   
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%% 2 Minkowski Function
    figure(12)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,1),'r');
    grid on; xlabel('r'); ylabel('La'); 
    title('Second Minkowski Function (q=0.001)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,1),'b');
    for i= 2:m

        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,3,i),'b');

    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,m+1),'r','LineWidth',1.5);   
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%% 3 Minkowski Function
    figure(13)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,1),'r');
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function (q=0.001)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,1),'b');
    for i= 2:m    
        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,4,i),'b');
    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,m+1),'r','LineWidth',1.5); 
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%%%% noise sensitivity
    figure(14)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('First Minkowski Function for a Noise Sensitivity of 0.001');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,m+1),'r','LineWidth',1.0);
    hold off

    figure(15)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Second Minkowski Function for a Noise Sensitivity of 0.001');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,m+1),'r','LineWidth',1.0);
    hold off

    figure(16)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function for a Noise Sensitivity of 0.001');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,m+1),'r','LineWidth',1.0);
    hold off
    fprintf('Task4a Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 4b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task4b()
    lambda=50; b=0.04; a1=0.2; a2=0.5; p=0.5; W=[0,1;0,1];
    m=20; k=40;
    spacing=0.005; q=0.01;  % Disturbance coefficient
    minkovFnTable = zeros(k+1,4,m+1);
    distMinkovFnTable = zeros(k+1,4,m+1);

    for i = 1:m    
        rzn = rBoolEllipse(lambda,b,a1,a2,p,W);  % Sample Realisation
        rbool =digitizeEllSys(rzn,W,spacing);  % digitization of the generated Realisation
        minkovFnTable(:,:,i) =estQMinkowskiFcts(rbool,k,spacing); 
        rbool1 = disturbImage(rbool,q); % Disturbing the image
        distMinkovFnTable(:,:,i) =estQMinkowskiFcts(rbool1,k,spacing);
    end

    %minkovFnTable(:,:)
    %Average of the Minkovski Functional
    for i =2:4
        for j = 1:k+1
            avg =0; avg1=0;
            for k =1:m
                temp = double( minkovFnTable(j,i,k));
                temp1 = double( distMinkovFnTable(j,i,k));
                avg = temp + avg;
                avg1 = temp1 + avg1;
               % minkovFnTable(j:i:21)= minkovFnTable(j:i:21)+ minkovFnTable(j:i:k);

            end
            minkovFnTable(j,i,m+1)= avg/m;
            distMinkovFnTable(j,i,m+1)= avg1/m;
        end
    end
    %%%%%%%%%%%%%%%%
    % file to create the plot for undisturbed minkovski fn table
    m = 20;

    figure(17)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,1),'r');
    grid on; xlabel('r'); ylabel('Aa'); 
    title('First Minkowski Function (q=0.01)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,1),'b');
    for i= 2:m

        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,2,i),'b');

    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,m+1),'r','LineWidth',1.5);   
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%% 2 Minkowski Function
    figure(18)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,1),'r');
    grid on; xlabel('r'); ylabel('La'); 
    title('Second Minkowski Function (q=0.01)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,1),'b');
    for i= 2:m

        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,3,i),'b');

    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,m+1),'r','LineWidth',1.5);   
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%% 3 Minkowski Function
    figure(19)
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,1),'r');
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function (q=0.01)');
    hold on
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,1),'b');
    for i= 2:m    
        plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,i),'r');
        plot(minkovFnTable(:,1,1),minkovFnTable(:,4,i),'b');
    end
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,m+1),'r','LineWidth',1.5); 
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'b','LineWidth',1.5);
    hold off
    %%%%%%%%% noise sensitivity
    figure(20)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,2,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('First Minkowski Function for a Noise Sensitivity of 0.01');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,2,m+1),'r','LineWidth',1.0);
    hold off

    figure(21)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,3,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Second Minkowski Function for a Noise Sensitivity of 0.01');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,3,m+1),'r','LineWidth',1.0);
    hold off

    figure(22)
    plot(minkovFnTable(:,1,1),minkovFnTable(:,4,m+1),'b','LineWidth',1.0);
    grid on; xlabel('r'); ylabel('Xa'); 
    title('Third Minkowski Function for a Noise Sensitivity of 0.01');
    hold on
    plot(distMinkovFnTable(:,1,1),distMinkovFnTable(:,4,m+1),'r','LineWidth',1.0);
    hold off
    fprintf('Task4b Completed ''\n')
    fprintf('Task4 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Envelope = task5()
    fprintf('Running Task5 ''\n')
    lambda=50; b=0.04; a1=0.2; a2=0.5; p=0.5; W=[0,1;0,1];
    m=20; k_max=40; spacing=0.005; 
    Envelope = calcEnvelope(lambda,a1,a2,b,p,W,spacing,m,k_max);
    fprintf('Task5 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% TASK 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function task6()
    fprintf('Running Task6 ''\n')
    %%%%%%%%%%%taska%%%%%%%%%%
    spacing=0.01;
    k=30;
    testStat = dlmread('data_Radhakrishnan.txt');
    minkovFnTable = estQMinkowskiFcts(testStat,k,spacing);

    %imshow(testStat);
    figure(23)
    plot(minkovFnTable(:,1),minkovFnTable(:,2),'b');
    grid on; xlabel('r'); ylabel('Aa');
    title('First Minkowski Function');

    figure(24)
    plot(minkovFnTable(:,1),minkovFnTable(:,3),'b');
    grid on; xlabel('r'); ylabel('La');
    title('Second Minkowski Function');

    figure(25)
    plot(minkovFnTable(:,1),minkovFnTable(:,4),'b');
    grid on; xlabel('r'); ylabel('Xa');
    title('Third Minkowski Function');
    fprintf('Task6a Completed ''\n')
    %%%%%%%%%%%taskb%%%%%%%%%%
    a1 = 0.3; a2 =1.2; b= 0.1;
    k= 30; spacing = 0.01;
    noOfPoints=15;

    %Reading from file 
    testStat = dlmread('data_Radhakrishnan.txt');
    MiTab = estQMinkowskiFcts(testStat,k,spacing);
    inGuess =[10;0.5];
    [modPar,fmin] = fminsearch(@getModelParam,inGuess,[],MiTab,noOfPoints, spacing,a1,a2,b);
    lambda = modPar(1,1);
    p = modPar(2,1);
    fprintf("Intensity factor = %f'\n'",lambda);
    fprintf("Probabilty(p) = %f'\n'",p);
    fprintf("fmin = %f'\n'",fmin);
    fprintf('Task6b Completed ''\n')
%%%%%%%%%%%taskc%%%%%%%%%%
    m = 39;                    %Number of plot
    W = [0,10;0,10];           %Window Size 

    ALXq = estQMinkowskiFcts(testStat,k,spacing);

    Aa_1(:,1) = ALXq(:,2);
    La_1(:,1) = ALXq(:,3);
    Xa_1(:,1) = ALXq(:,4);

    Envelope = calcEnvelope(lambda,a1,a2,p,b,W,spacing,m,k);
    
    L1 = Envelope(:,1);
    U1 = Envelope(:,2);
    L2 = Envelope(:,3);
    U2 = Envelope(:,4);
    L3 = Envelope(:,5);
    U3 = Envelope(:,6);
    r = spacing.*(0:k)';

    %Plotting "r V/s A(r)"
    figure(26);
    hold on;
    title('r v/s A(r)')
    xlabel('r')
    ylabel('A(r)')
    plot(r,U1,'b.',r,Aa_1,'r-',r,L1,'g.')
    hold off;

    %Plotting "r V/s L(r)"
    figure(27);
    hold on;
    title('r v/s L(r)')
    xlabel('r')
    ylabel('L(r)')
    plot(r,U2,'b.',r,La_1,'r-',r,L2,'g.')
    hold off;

    %Plotting "r V/s X(r)"
    figure(28);
    hold on;
    title('r v/s X(r)')
    xlabel('r')
    ylabel('X(r)')
    plot(r,U3,'b.',r,Xa_1,'r-',r,L3,'g.')
    hold off;
    fprintf('Task6c Completed ''\n')

%%%%%%%%%%%taskd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = 99;                     %Number of model realizations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deviance between values from realisations and expected value,
    % Data_Radhakrishnan.txt and expected value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:m
        
        [Matrix] = rBoolEllipse(lambda,a1,a2,p,b,W);
        [B] = digitizeEllSys(Matrix,W,spacing);
        ALXq = estQMinkowskiFcts(B,k,spacing);
        % Only one Minkowski function is required
        A(:,i) = ALXq(:,2);
        
    end
    
    %Storing the values of Minkowski functions in a new variable
    A_realisations = A;
    
    %Calculation of the Expected Function
    r = spacing.*(0:k)';
        
    for i=1:k+1
        A_expected(i,1) = 1-exp(-lambda*(r(i)*r(i) + (2*r(i)/pi)*((p*function_L(a1,b)) + (1-p)*function_L(a2,b))+ p*pi*a1*b + (1-p)*pi*a2*b));
    end
    
    % Deviance
    delt_realisations = zeros(size(A_realisations));
    for i = 1:size(A_realisations,2)
        delt_realisations(:,i) = (A_realisations(:,i) - A_expected).^2;
    end
    
    DeltaT_realisations = sum(delt_realisations);
    Store_DeltaT(:,1) = sort(DeltaT_realisations);
    
    testStat = dlmread('data_Radhakrishnan.txt');
    ALXq = estQMinkowskiFcts(testStat,k,spacing);
    
    %Copying the values from ALXq
    A_datapoints(:,1) = ALXq(:,2);
    
    % Deviance
    delt_datapoints = zeros(size(A_datapoints));
    for i = 1:size(A_datapoints,2)
        delt_datapoints(:,i) = (A_datapoints(:,i) - A_expected).^2;
    end
    
    DeltaT_datapoints = sum(delt_datapoints);
    Store_DeltaT(:,2) = sort(DeltaT_datapoints);
    
    
    dlmwrite('Task6D_expected.txt', Store_DeltaT);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Deviance between values from realisations and mean values,
    % Data_Radhakrishnan.txt and mean values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clc;
    clear all;
    close all;
    W = [0 10; 0 10];           %Window size
    spacing = 0.01;             %Side length of pixels
    k = 20;                     %Maximum range (measured in numbers of pixels)
    m = 99;                     %Number of model realizations
    
    %Values from task 6(b)
    lambda = 6.9401;            
    p = 0.4185;
    a1=0.3;
    a2=1.2;
    b=0.1;
    
    A_mean = zeros(k+1,1);
    
    for i=1:m
    
        [Matrix] = rBoolEllipse(lambda,a1,a2,p,b,W);
        [B] = digitizeEllSys(Matrix,W,spacing);
        ALXq = estQMinkowskiFcts(B,k,spacing);
        
        % Only one Minkowski function is required
        A(:,i) = ALXq(:,2);
        
    end
    
    % Mean value of the Minkowski function
    for i=1:k+1
        A_sum = 0;
        for j=1:m
            A_sum = A_sum + A(i,j);
        end
        A_mean(i,1) = A_sum/m;
    end 
    
    %Storing the values of Minkowski functions
    A_realisations = A;
    
    % Deviance
    delt_realisations = zeros(size(A_realisations));
    for i = 1:size(A_realisations,2)
        delt_realisations(:,i) = (A_realisations(:,i) - A_mean).^2;
    end
    
    DeltaT_realisations = sum(delt_realisations);
    Store_DeltaT(:,1) = sort(DeltaT_realisations);
    
    
    testStat = dlmread('data_Radhakrishnan.txt');
    ALXq = estQMinkowskiFcts(testStat,k,spacing);
    
    %Copying the values from ALXq
    A_datapoints(:,1) = ALXq(:,2);
    
    % Deviance
    delt_datapoints = zeros(size(A_datapoints));
    for i = 1:size(A_datapoints,2)
        delt_datapoints(:,i) = (A_datapoints(:,i) - A_mean).^2;
    end
    
    DeltaT_datapoints = sum(delt_datapoints);
    Store_DeltaT(:,2) = sort(DeltaT_datapoints);
    
    dlmwrite('Task6D_mean.txt', Store_DeltaT);

    fprintf('Task6d Completed ''\n')
    fprintf('Task6 Completed ''\n')
end

%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%