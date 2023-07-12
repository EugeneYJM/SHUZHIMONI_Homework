clear all;close all;clc;
L=1; % meter
Nx=100; %number of grid point 
x=linspace(0,L,Nx); %x grid
dx=L/(Nx-1); %x grid resolution, 1 cm 
K=1.e-2; %thermal diffuction coefficient , m^2/s
A=5;  %initial temp. amp. C degree.
dt=0.5*dx^2/K; % time step required by the CFD stability condition:
mu= dt*K/dx^2 ;%<=0.5
lambda=pi/L;  %initial wave number 
nT= 10000; %total number of iteration
%%Exact solution  normal distribution
figure; 
T0=A*exp(-(x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
an=2/Nx*fft(T0)';
kn=[0:Nx/2  , -Nx/2+1:-1]'.*2*pi/L; %L=Nx*dx
expikx= exp(1j.*(kn*x));
pic_num=1;
for it=0:200:nT
    T_exact = (an.*exp(-K*kn.^2*it*dt))'*expikx ;
    area(x, real(T_exact));
    title(['Exact solution ',num2str(it*dt),' s']);
    xlim([0 L]);
    ylim([0  A]);
    pause(0.1)
        % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if pic_num==1
        imwrite(imind,cm,'assignment2_1d_diff_FTCS_Exact solution_normaldistribution.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'assignment2_1d_diff_FTCS_Exact solution_normaldistribution.gif','gif','WriteMode','append');
    end
    pic_num=pic_num+1;
end