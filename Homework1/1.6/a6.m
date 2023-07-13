%% assignment 6.1  对流扩散 显式FTCS
clear all;close all;clc;
L=10; % meter
Nx=1000; %number of grid point 
x=linspace(0,L,Nx); %x grid
dx=L/(Nx-1); %x grid resolution, 1 cm 
A=100;  %initial temp. amp. C degree.
nT= 200; %total number of iteration
u=1;%m/s
K=1.e-2; %thermal diffuction coefficient , m^2/s
dt=0.5*dx^2/K; % time step required by the CFD stability condition:
mu= (dt*K)/(dx^2) ;%<=0.5
lambda=(u*dt)/dx;%%0.5  courant<1

T0=A*exp(  - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); %t
xi=2:Nx-1;

figure; 
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<2000
    count=count+1;
    T1(xi)=(1-2*mu)*T0(xi)+(mu-lambda/2)*T0(xi+1)+(mu+lambda/2)*T0(xi-1);
    T1(1)=T1(2);%boundary condition：dT/dx(x=0,x=L)=0
    T1(end)=T1(end-1);%boundary condition
    T0=T1;
    if mod(count,25)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['FTCS 1D convection eq. ',num2str(count*dt),' s',' Lambda: ',num2str(lambda)]);
        xlabel('x')
        ylabel('T')
        pause(0.025);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 25
            imwrite(imind,cm,'assignment6.1_1d_conv_FTCS.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment6.1_1d_conv_FTCS.gif','gif','WriteMode','append');
        end
    end
end


%% assignment 6.2  对流扩散 隐式BTCS
clear all;close all;clc;
L=10; % meter
Nx=1000; %number of grid point
x=linspace(0,L,Nx); %x grid
dx=L/(Nx-1); %x grid resolution, 1 cm
A=100;  %initial temp. amp. C degree.
nT= 500; %total number of iteration
u=1;%m/s
K=1.e-2; %thermal diffuction coefficient , m^2/s
dt=2*dx^2/K; % time step required by the CFD stability condition:
mu= dt*K/dx^2 ;%<=0.5
lambda=u*dt/dx;%% courant>1

T0=A*exp( - (x-0.5*L).^2/(0.1*L)^2  );%%initial condition__normal distribution
T1=zeros(size(T0)); %t
xi=2:Nx-1;

D=zeros(Nx-2,Nx-2);
for i=1:Nx-2
    D(i,i)= - (1+2*mu);
    if i<Nx-2
        D(i,i+1)=mu-lambda/2;
    end
    if i>1
        D(i,i-1)=mu+lambda/2;
    end
end
D(1,1)=-1-mu+lambda/2;
D(end,end)=-1-mu-lambda/2;

D_inv= D^-1;
figure;
drawnow;
set(gca,'xlim',[0 L],'ylim',[0 A])
count=0;
while count<nT
    count=count+1;
    T1(xi)= - D_inv*T0(xi)';
    T1(1)=T1(2);
    T1(end)=T1(end-1);
    T0=T1;
    if mod(count,10)==0
        area(x,T1,-A);
        xlim([0 L]);
        ylim([-A A]);
        title(['BTCS 1D convection eq. ',num2str(count*dt),' s',' Lambda: ',num2str(lambda)]);
        xlabel('x')
        ylabel('T')
        pause(0.025);
        % Capture the plot as an image
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if count == 10
            imwrite(imind,cm,'assignment6.2_1d_conv_BTCS.gif','gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,'assignment6.2_1d_conv_BTCS.gif','gif','WriteMode','append');
        end
    end
end

