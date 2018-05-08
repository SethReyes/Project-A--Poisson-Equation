%%Poisson Equation







% Domain of Interest: ax < x < bx, ay < y < by
ax = -pi; bx = pi;
ay = -pi; by = pi;
% Grid
N=100; step = 1/N;
xLen=2*pi; yLen=2*pi; % Square domain


x=linspace(ax,bx,N);
y=linspace(ay,by,N);
[xx,yy]=meshgrid(x,-y);
h=(bx-ax)/(N-1);

%% Boundary Conditions
U = zeros(N,N);
fb = (by-yy).^2*cos(pi.*yy/by);
gb = (by-yy).^2.*yy;
% Walls of U 
U(1:N-1,1) = fb(1:N-1,1);
U(1:N-1,N) = gb(1:N-1,N);
temp = (by-ay)^2*cos(pi*ay/by) + (xx-ax)/(bx-ax)*((by-ay)^2*ay-(by-ay)^2*cos(pi*ay/by));
U(N,1:N) = temp(N,1:N);

F = zeros(N,N);
F = cos(pi/2*(2.*((xx-ax)/(bx-ax))+1)).*sin(pi.*(yy-ay)/(by-ay));

%% Gauss-Siedel Method
GS_U=U;
maxError=1;
GSiterations=0;
tic;
while 1e-7 < maxError
    prevU=GS_U;
    % Neumann Boundary Condition 2*Udown+Uleft+Uright
    GS_U(1,2:N-1) = 1/4*(2*GS_U(1+1,2:N-1)+GS_U(1,(2:N-1)-1)+GS_U(1,(2:N-1)+1)+h^2*F(1,2:N-1));
    for i=2:N-1
        for j=2:N-1  
            GS_U(i,j) = (1/4)*(prevU(i+1,j)+GS_U(i-1,j)+prevU(i,j+1)+GS_U(i,j-1)+(h^2)*F(i,j));
        end
    end
    %Error check
    maxError=0;
    for j=1:N-1
        for i=2:N-1
            difference=abs(prevU(i,j)-GS_U(i,j));
            if difference > maxError
                maxError=difference;
            end
        end
    end
    GSiterations=GSiterations+1;
end
[Uxx,Uyy]=gradient(-GS_U,h,h);
tU = sqrt(Uxx.^2+Uyy.^2);

GSiterations
GS_U;

% Vector field of U
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
%Countour of U
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

GS_time=toc
%% Successive Over Relaxation
tic;
SOR_U=U;
maxError=1;
SOR_iterations=0;
w=1.979
while 1e-7 < maxError
    prevU=SOR_U;
    % Neumann Boundary Condition 2*Udown+Uleft+Uright
    SOR_U(1,2:N-1) = 1/4*(2*SOR_U(1+1,2:N-1)+SOR_U(1,(2:N-1)-1)+SOR_U(1,(2:N-1)+1)+h^2*F(1,2:N-1));
    for i=2:N-1
        for j=2:N-1  
        SOR_U(i,j) = (w/4)*(prevU(i+1,j)+SOR_U(i-1,j)+prevU(i,j+1)+SOR_U(i,j-1)+(h^2)*F(i,j))+(1-w)*prevU(i,j);
        end
    end
    %Error check
    maxError=0;
    for j=1:N-1
        for i=2:N-1
            difference=abs(prevU(i,j)-SOR_U(i,j));
            if difference > maxError
                maxError=difference;
            end
        end
    end
    SOR_iterations=SOR_iterations+1;
end
[Uxx,Uyy]=gradient(-SOR_U,h,h);
tU = sqrt(Uxx.^2+Uyy.^2);

SOR_iterations
SOR_U;

% Vector field of U
figure(3)
set(gcf,'units','normalized','position',[0.2 0.1 0.3 0.32]);
surf(xx,yy,SOR_U);
xlabel('x'); ylabel('y'); zlabel('U');
title('Solution for U using the SOR Method','fontweight','normal');
rotate3d
box on
axis tight
k =  colorbar;
k.Label.String = 'U';
colormap default;
%Countour of U
figure(4)
set(gcf,'units','normalized','position',[0.5 0.1 0.3 0.32]);
contourf(xx,yy,SOR_U,25);
xlabel('x '); ylabel('y ');
shading interp
title('Contour of U using the SOR Method','fontweight','normal');
box on
k =  colorbar;
k.Label.String = 'U';
axis square
box on

SOR_time=toc

