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
U(2:N-1,1) = fb(2:N-1,1);
U(2:N-1,N) = gb(2:N-1,N);
temp = (by-ay)^2*cos(pi*ay/by) + (xx-ax)/(bx-ax)*((by-ay)^2*ay-(by-ay)^2*cos(pi*ay/by));
U(N,2:N-1) = temp(N,2:N-1);
% Corners of U
U(N,1) = (U(N-1,1)+U(N,2))/2;
U(N,N) = (U(N,N-1)+U(N-1,N))/2;
F = zeros(N,N);
F = cos(pi/2*(2.*((xx-ax)/(bx-ax))+1)).*sin(pi.*(yy-ay)/(by-ay));


%% Guassian Siedel
maxError=1;
GaussSiedeliterations=0;
while 1e-7 < maxError
    prevU=U;
    U(1,1) = (U(1,2) + U(2,1))/2;
    U(1,N) = (U(1,N-1)+U(2,N))/2;
    % Neumann Boundary Condition 2*Udown+Uleft+Uright
    U(1,2:N-1) = 1/4*(2*U(1+1,2:N-1)+U(1,(2:N-1)-1)+U(1,(2:N-1)+1)+h^2*F(1,2:N-1));
    for i=2:N-1
        for j=2:N-1  
        U(i,j) = (1/4)*(prevU(i+1,j)+prevU(i-1,j)+prevU(i,j+1)+prevU(i,j-1)+(h^2)*F(i,j));
        end
    end
    %Error check
    maxError=0;
    for i=1:N-1
        for j=2:N-1
            difference=abs(prevU(i,j)-U(i,j));
            if difference > maxError
                maxError=difference;
            end
        end
    end
    GaussSiedeliterations=GaussSiedeliterations+1;
end
[Uxx,Uyy]=gradient(-U,h,h);
tU = sqrt(Uxx.^2+Uyy.^2);

GaussSiedeliterations
U;

% Vector field of U
figure(1)
mesh(xx,yy,U);
xlabel('x'); ylabel('y'); zlabel('U');
title('Solution using Gauss Seidel Method','fontweight','normal');
rotate3d
box on
axis tight
h =  colorbar;
h.Label.String = 'U ';
