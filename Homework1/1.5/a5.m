clear all;close all;clc;
L=1; % meter
Nx=100; %number of grid point 
x=linspace(0,L,Nx); %x grid
dx=L/(Nx-1); %x grid resolution, 1 cm 
A=5;  %initial temp. amp. C degree.
dt=0.8*dx; % time step required by the CFD stability condition:
nT= 10000; %total number of iteration
u=1;%m/s
lambda=u*dt/dx;%%0.8  courant <1

%-----------------------------------------------------------------------------------------------------
%%%%%%%蛙跃
%FTCS: Euler forward in time central in space:
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); 
T2=zeros(size(T0)); 
xi=2:Nx-1;

T1(xi)=T0(xi)-lambda/2*(T0(xi+1)-T0(xi-1));%FTCS
T1(1)=T1(end-1);%boundary condition：T(1,t)=T(n-1,t)
T1(end)=T1(2);%boundary condition:T(n,t)=T(2,t)

figure; %numerical solution by CTCS scheme
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T2(xi)=T0(xi)-lambda*(T1(xi+1)-T1(xi-1));%%ctcs
    T2(1)=T2(end-1);%boundary condition：T(1,t)=T(n-1,t)
    T2(end)=T2(2);%boundary condition:T(n,t)=T(2,t)
    
    T0=T1; %update T0, T1 to predict T2 for the next time step
    T1=T2;  %
    T1(1)=T1(end-1);%boundary condition：T(1,t)=T(n-1,t)
    T1(end)=T1(2);%boundary condition:T(n,t)=T(2,t)

    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['CTCS 1D convection eq',num2str(count*dt),' s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment5_1d_conv_CTCS_bigdt.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment5_1d_conv_CTCS_bigdt.gif','gif','WriteMode','append');
        end
    end 
end

%------------------------------------------------------------------------------
%%%%%%%迎风（u>0）
T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); %t
xi=2:Nx-1;

figure; 
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T1(xi)=T0(xi)-lambda*(T0(xi)-T0(xi-1));
    T1(1)=T1(end-1);%boundary condition：T(1,t)=T(n-1,t)
    T1(end)=T1(2);%boundary condition:T(n,t)=T(2,t)
    T0=T1;
    if mod(count,200)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['FTFS 1D convection eq. ',num2str(count*dt),' s']);
        xlabel('x')
        ylabel('T')
        pause(0.1);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 200
            imwrite(imind,cm,'assignment5_1d_conv_FTFS_bigdt.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment5_1d_conv_FTFS_bigdt.gif','gif','WriteMode','append');
        end
    end
end
