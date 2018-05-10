%%  Project A - Poisson Equation
%   Scientific Computing for Mechanical Engineers
%   MECE 5397
%   Version: APc1-3
%   Seth Reyes
%   5/9/2018
 

%   This code uses two methods to find the solution of the Poisson Equation.
%   The Gauss-Seidel method and Successive Over Relaxation method is used
%   and are compared in the end. The code prompts the user for the number 
%   of steps points the user wants the program to solve. This number 'N'
%   will be used for both the X direction and the Y direction to form a 2-D 
%   grid to solve. Next, the user will be prompted to input a value for
%   the error which will be in the order of 10's.
clc, clear all, close all


%% Prompts and Reboot Check
GScompleted=0;  % Will run GS computations if set to 0
reboot='N';     % assumes no reboot
GS_time=0;     % If rebooted, storeGS_time will be added onto tic toc
SOR_time=0; %Initial time = 0
if exist('checkpoint.mat', 'file') == 2     %checks for checkpoint from last run
    acceptableResponse=0;   %acceptable reponse initialized as false
    disp('It has been detected that the previous run has crashed or has been terminated early.')
    while acceptableResponse == 0
        reboot = input('Would you like to start where you left off? (Y/N)\n\n','s'); % self explanatory
        if isequal(reboot,'Y') || isequal(reboot,'y')
            load checkpoint.mat %loads in last checkpoint
            disp(' ')
            disp('Loading...')
            acceptableResponse=1;   % Response Valid
        elseif reboot == 'N' || reboot == 'n'
            acceptableResponse=1;   % Response Valid
            disp(' ')
        else
            disp('Input Error') % Response Invalid
        end
    end    
end
if reboot == 'N' || reboot == 'n'   % If no reboot is required.. Prompts. Else start at calculations
    disp('Select the value of F to solve the Poisson Equation for.')
    disp(' ')
    disp('Case #1: F = cos(pi/2*(2.*((xx-ax)/(bx-ax))+1)).*sin(pi.*(yy-ay)/(by-ay))')
    disp('Case #2: F = 0')
    acceptableResponse=0;
    while acceptableResponse == 0 %checks for an acceptable reponse to the question
        solveF=input('\nSelecting Case #');
        if solveF == 1 || solveF == 2
            acceptableResponse=1;   % Response Valid
        else
            disp('Input Error') % Response Invalid
        end
    end
    acceptableResponse=0;
    while acceptableResponse == 0
        watch=input('\n\nWould you like to watch the grid convergence? (Y/N)\nWarning: Selecting "Y" can severely dimiminish performance. Please\n         choose numbers smaller than the recommended values.\n\n','s');
        if watch == 'Y' || watch == 'y' || watch == 'N' || watch == 'n'
            acceptableResponse=1;% Response Valid
        else
            disp('Input Error') % Response Invalid
        end
    end
    disp(' ')
    disp(' ')
    disp('Select the number of steps that you would like to solve for. "N" will be used in the X and Y direction.')
    disp('"N" must be an integer greater than 2, but larger numbers can be quite taxing on the computer.')
    disp('Recommended: N = 100')
    N=input('\nValue of N = ');
    disp(' ')
    disp('Please choose the level of accuracy (10^-X) to solve for. Larger values of "X" can be quite taxing on the computer.')
    disp('Recommended: X = 7')
    acceptableResponse=0;
    while acceptableResponse == 0 %checks for an acceptable reponse to the question    
        Err=input('\nValue of X = ');
        if isnumeric(Err)
            acceptableResponse=1;   % Response Valid
        else
            disp('Input Error')% Response Invalid
        end
    end
    fprintf('Accuracy = 10^-%1.0f\n\n(To terminate run, press Ctrl + C)\n\nLoading...\n\n',Err)

    %% Boundary Conditions
    % Domain of Interest: ax < x < bx, ay < y < by
    ax = -pi; bx = pi;
    ay = -pi; by = pi;
    % Grid
    xLen=2*pi; yLen=2*pi; % Square domain
    x=linspace(ax,bx,N);    % Creates array from ax to bx
    y=linspace(ay,by,N);    % Creates array from ay to by
    [xx,yy]=meshgrid(x,-y); % Creates matrix of above arrays of only that vector
    h=(bx-ax)/(N-1);    % Step size
    U = zeros(N,N);     % Empty U matrix
    fb = (by-yy).^2*cos(pi.*yy/by); %solves fb to be put into the B.C.
    gb = (by-yy).^2.*yy;            %solves gb to be put into the B.C.
    % Boundaries of U 
    U(1:N-1,1) = fb(1:N-1,1); %Left Boundary Condition
    U(1:N-1,N) = gb(1:N-1,N); %Right Boundary Condition
    temp = (by-ay)^2*cos(pi*ay/by) + (xx-ax)/(bx-ax)*((by-ay)^2*ay-(by-ay)^2*cos(pi*ay/by)); %temp matrix for B.C.
    U(N,1:N) = temp(N,1:N); %Bottom B.C.
    F = zeros(N,N); %For F=0
    if solveF==1
        F = cos(pi/2*(2.*((xx-ax)/(bx-ax))+1)).*sin(pi.*(yy-ay)/(by-ay)); % For F = prescribed condition
    end
%% Gauss-Siedel Method
    GS_U=U; %Stores U as GS_U without overwriting U
    maxError=1; % temp initial error
    GS_iterations=0; % start counting iterations
end % End of "No Reboot" stage

storeGS_time=GS_time;   % If rebooted, storeGS_time will be added onto tic toc
if GScompleted==0 %If Gausss-Seidel completed last run, then skip G-S Calculations
    tic;        %time start
    while 10^-(Err) < maxError% while error is bigger than prescribed accuracy
        prevU=GS_U; % Saving Previous U
        for j=2:N-1
            for i=2:N-1  
                GS_U(i,j) = (1/4)*(prevU(i+1,j)+GS_U(i-1,j)+prevU(i,j+1)+GS_U(i,j-1)+(h^2)*F(i,j)); %GS_U
            end
            % Neumann Boundary Condition 2*Udown+Uleft+Uright
            GS_U(1,j) = 1/4*(2*GS_U(2,j)+GS_U(1,j-1)+prevU(1,j+1)+h^2*F(1,j));            
        end
        %Error check
        maxError=max(max(abs(prevU-GS_U))); % finds the error that is greatest between nodes of previous and new
        GS_time=toc;    % temporarily saves toc for checkpointing
        GS_iterations=GS_iterations+1; %iterations count 1 up
        if mod(GS_iterations,1000)==0   % Checkpoint: every 1000 iterations, save
            save checkpoint.mat N Err GS_U prevU GS_iterations GS_time maxError h F U xx yy %Save variables...
        end
        % Animation for Grid Convergence
        if watch=='y' || watch=='Y'
            figure(1)
            set(gcf,'units','normalized','position',[0.2 0.5 0.3 0.32]);        
            surf(xx,yy,GS_U);         
            xlabel('x'); ylabel('y'); zlabel('U');
            title('Solution for U using the Gauss-Seidel Method','fontweight','normal');
            rotate3d
            box on
            axis tight
            k =  colorbar;
            k.Label.String = 'U';
            colormap default;
        end
    end
GS_time=toc+storeGS_time; %adds previous run time. If not reboot, then storeGS_time=0
end
GScompleted=1; % If reboot, won't calculate GS again

%% Successive Over Relaxation

storeSOR_time=SOR_time; %if no checkpoint for SOR, = 0, else storing previous run
tic; %start timer for SOR
SOR_U=U; %Store U for SOR without overwriting U
maxError=1; % temp false error
SOR_iterations=0; % SOR iterations
w=round(2/(1+sin(h/2)),2); % Initial guess for w. Better w correspond with different N values
while 10^-(Err) < maxError
    prevU=SOR_U;    %Stores Previous U from last calculated U    
    for j=2:N-1
        for i=2:N-1  
            SOR_U(i,j) = (w/4)*(prevU(i+1,j)+SOR_U(i-1,j)+prevU(i,j+1)+SOR_U(i,j-1)+(h^2)*F(i,j))+(1-w)*prevU(i,j); %solving for U(i,j) using SOR
        end
        % Neumann Boundary Condition 2*Udown+Uleft+Uright
        SOR_U(1,j) = 1/4*(2*SOR_U(2,j)+SOR_U(1,j-1)+prevU(1,j+1)+h^2*F(1,j));        
    end
    %Error check
    maxError=max(max(abs(prevU-SOR_U)));    % finds the error that is greatest between nodes of previous and new
    SOR_iterations=SOR_iterations+1;        % iterations count 1 up
    SOR_time=toc;                           % temporarily saves toc for checkpointing
    if mod(SOR_iterations,1000)==0          % Checkpoint: every 1000 iterations, save
        save checkpoint.mat N Err GS_U GS_iterations GS_time SOR_U prevU SOR_iterations GScompleted maxError h F U xx yy SOR_time %saves variables
    end
    % Animation for Grid Convergence
    if watch == 'y' || watch == 'Y'
        figure(3)
        set(gcf,'units','normalized','position',[0.2 0.1 0.3 0.32]);
        surf(xx,yy,SOR_U);
        xlabel('x'); ylabel('y'); zlabel('U');
        title('Solution for U using the Successive Over Relaxation Method','fontweight','normal');
        rotate3d
        box on
        axis tight
        k =  colorbar;
        k.Label.String = 'U';
        colormap default;
    end
end

SOR_time=toc+storeSOR_time; % Adds previous run time to current run time

%% Plots and Graphs

% Vector field of U using GS
figure(1)
set(gcf,'units','normalized','position',[0.2 0.5 0.3 0.32]);
if N>99
    mesh(xx,yy,GS_U);
else
    surf(xx,yy,GS_U); % surf gets ugly for N>100
end
xlabel('x'); ylabel('y'); zlabel('U');
title('Solution for U using the Gauss-Seidel Method','fontweight','normal');
rotate3d
box on
axis tight
k =  colorbar;
k.Label.String = 'U';
colormap default;

%Countour of U using GS
figure(2)
set(gcf,'units','normalized','position',[0.5 0.5 0.3 0.32]);
contourf(xx,yy,GS_U,25);
xlabel('x '); ylabel('y ');
shading interp
title('Contour of U using the Gauss-Seidel Method','fontweight','normal');
box on
k =  colorbar;
k.Label.String = 'U';
axis square
box on

% Vector field of U using SOR
figure(3)
set(gcf,'units','normalized','position',[0.2 0.1 0.3 0.32]);
if N>99
    mesh(xx,yy,SOR_U);
else
    surf(xx,yy,SOR_U);
end
xlabel('x'); ylabel('y'); zlabel('U');
title('Solution for U using the Successive Over Relaxation Method','fontweight','normal');
rotate3d
box on
axis tight
k =  colorbar;
k.Label.String = 'U';
colormap default;

%Countour of U using SOR
figure(4)
set(gcf,'units','normalized','position',[0.5 0.1 0.3 0.32]);
contourf(xx,yy,SOR_U,25);
xlabel('x '); ylabel('y ');
shading interp
title('Contour of U using the Successive Over Relaxation Method','fontweight','normal');
box on
k =  colorbar;
k.Label.String = 'U';
axis square
box on

%% Results

fprintf('\n\nThe amount of time it took to find a solution for U using the Gauss-Seidel Method is: %5.4f seconds',GS_time)
fprintf('\nThe amount of time it took to find a solution for U using the Successive Over Relaxation Method is: %5.4f seconds',SOR_time)
fprintf('\nThe number of iterations it took to find a solution for U using the Gauss-Seidel Method is: %5.0f',GS_iterations)
fprintf('\nThe number of iterations it took to find a solution for U using the Successive Over Relaxation Method is: %5.0f',SOR_iterations)
fprintf('\nThe total sum of the difference between the %1.0f nodes of the G-S Method and the SOR Method is: %5.4f\n\n\n',(N^2),sum(sum(abs(GS_U-SOR_U))));

delete checkpoint.mat %removes checkpoints if code runs all the way through