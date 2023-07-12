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
%%Exact solution
figure;%
drawnow;
pic_num=1;
for it=0:200:nT
    T_exact = A*exp(-K*lambda^2*it*dt)*cos(lambda*x);
    area(x,T_exact, -A);
    title(['Exact solution ',num2str(it*dt),' s']);
    xlim([0 L]);
    ylim([-A A]);
    pause(0.1)
    % Capture the plot as an image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if pic_num==1
        imwrite(imind,cm,'assignment1_1d_diff_FTCS_Exact solution.gif','gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,'assignment1_1d_diff_FTCS_Exact solution.gif','gif','WriteMode','append');
    end
    pic_num=pic_num+1;
end

%----------------------------------------------------------------------------------------------------------------
%numerical solution   
T0=A*cos(1*x/L*pi); %%initial condition
T1=zeros(size(T0)); 
xi=2:Nx-1;

figure; %numerical solution by FTCS scheme:Euler forward in time central in space:
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+ mu*(  T0(xi+1)+T0(xi-1)  );%FTCS
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1; 
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['Numerical solution ',num2str(count*dt),' s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment1_1d_diff_FTCS_Numerical solution.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment1_1d_diff_FTCS_Numerical solution.gif','gif','WriteMode','append');
        end
    end 
end
