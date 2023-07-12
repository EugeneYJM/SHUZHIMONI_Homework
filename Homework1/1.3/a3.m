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
%% assignment 3.1  显性FTCS和隐性BTCS
%FTCS
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); %t
xi=2:Nx-1;

figure; %numerical solution by FTCS scheme
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+ mu*(T0(xi+1)+T0(xi-1)  );
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1;
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['FTCS 1D diffusion eq. ',num2str(count*dt),' s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment3.1_1d_diff_FTCS_explicit.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment3.1_1d_diff_FTCS_explicit.gif','gif','WriteMode','append');
        end
    end
end

%BTCS
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );
T1=zeros(size(T0)); %create an empty matrix for n+1 time layer
xi=2:Nx-1;  % innner grid points
% construct the differencing matrix
% BTCS: Euler Backword in time central differencing in space:
D=zeros(Nx-2,Nx-2);
for i=1:Nx-2
    D(i,i)= - (1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu;
    end
    if i>1
        D(i,i-1)=mu;
    end
end
D(1,1)=-(1+mu);
D(end,end)=-(1+mu);% T_x=0 in the implicit Martrix form
D_inv= D^-1;
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<10000
    count=count+1;
    T1(xi)= - D_inv*T0(xi)';
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1;
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['BTCS 1D diffusion eq. ' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment3.1_1d_diff_BTCS_implicit.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment3.1_1d_diff_BTCS_implicit.gif','gif','WriteMode','append');
        end
    end
end
%% assignment3.2 显性FTCS和隐性BTCS  时间步长dt增大
dt=0.8*dx^2/K; % time step required by the CFD stability condition:
%FTCS
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); %t
xi=2:Nx-1;

figure; %numerical solution by FTCS scheme
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+ mu*(T0(xi+1)+T0(xi-1)  );
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1;
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['FTCS 1D diffusion eq 1. ',num2str(count*dt),' s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment3.2_1d_diff_FTCS_explicit_1.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment3.2_1d_diff_FTCS_explicit_1.gif','gif','WriteMode','append');
        end
    end
end

%BTCS
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );
T1=zeros(size(T0)); %create an empty matrix for n+1 time layer
xi=2:Nx-1;  % innner grid points
% construct the differencing matrix
% BTCS: Euler Backword in time central differencing in space:
D=zeros(Nx-2,Nx-2);
for i=1:Nx-2
    D(i,i)= - (1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu;
    end
    if i>1
        D(i,i-1)=mu;
    end
end
D(1,1)=-(1+mu);
D(end,end)=-(1+mu);% T_x=0 in the implicit Martrix form
D_inv= D^-1;
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<10000
    count=count+1;
    T1(xi)= - D_inv*T0(xi)';
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1;
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['BTCS 1D diffusion eq 1. ' num2str(count*dt),'s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment3.2_1d_diff_BTCS_implicit_1.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment3.2_1d_diff_BTCS_implicit_1.gif','gif','WriteMode','append');
        end
    end
end