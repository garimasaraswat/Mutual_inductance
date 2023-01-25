% Program to find the theoretical value of Mutual inductance: calculated by
% using both coil dimensions, distance beween film and coils, Drive current
% and partitioning the film into N concentric annular rings and trial value
% of lambda ::::::::::::: GARIMA SARASWAT, 2015::::::::::::::::::

%----------------------------Initialising variables:'Please enter all values in S.I. unit
clear
%msgbox(,'UNITS,'warn');
ndh=18;%input('turns in each clockwise/anticlockwise section of drive coil');
npr=4;%input('layers in pickup coil');
nph=41;%input('turns in each layer of pickup coil');
R=input('Radius of film');
d=input('thickness of film');
sb=input('thickness of substrate');
wire=5e-05;%input('distance between each loop');
r_drive=0.9825e-3;%input('radius of drive coil');
r_pickup=0.9825e-3;%input('radius of pickup coil');
h_drive=0.61e-3;%input('seperation of nearest turn of drive coil to the film');
h_pickup=0.53e-3;%input('seperation of nearest turn of pickup coil to the film');
N=input('number of partitioned rings in the film');
%freq=input('Frequency of drive current');
Id=input('amplitude of drive current');
lambdain=input('trial value of lambda');
mu_o=4*pi*1e-7; 

% %---------------------------------------------------------------------------------------------------%
% %1. Dividing sample film into rings (using equal width) & effective thickness
% %---------------------------------------------------------------------------------------------------%
ls=1;
delta_r=R/N;
for i=1:N
    r_ring(i)=((i-1)*delta_r)+(delta_r/2);
end
for lambda=lambdain:10e-10:lambdain+100e-10
%for lambda=lambdain:1e-10:lambdain+100e_compl
    clear M k1square ksquare K E z1 K1 K2 E1 E2 z1 z2 B B1 B2 Ap_d1 Ap_d2 Ap_d Ap_r d_eff C A B J
    d_eff=lambda*sinh(d/lambda);
% %----------------------------------------------------------------------------------------------------%    
% %----------------------2. b_i = vector potential at i_th ring of film due to drive coil
% %----------------------------------------------------------------------------------------------------%
B=zeros(N,1); %intialising column matrix
B1=zeros(N,1);%for clockwise turns
B2=zeros(N,1);% for anti clockwise turns
for i=1:N
    for j=1:ndh %turns in each clockwise/anticlockwise section of drive coil 
        clear k1square z1 z2 K1 E1 K2 E2;
        z1=sb+h_drive+((j-1)*wire);%distance of clockwise loop from ring
        z2=sb+h_drive+((j+17)*wire);%distance of anti-clockwise loop from ring
    k1square=(4*r_drive*r_ring(i))/((r_drive+r_ring(i))^2 +z1^2);
    k2square=(4*r_drive*r_ring(i))/((r_drive+r_ring(i))^2 +z2^2);
    [K1,E1]=ellipke(k1square);
    [K2,E2]=ellipke(k2square);
    B1(i)=B1(i)+ ((1/pi)*r_drive*((2-k1square)*K1 - 2*E1)/(k1square*(sqrt((r_drive+r_ring(i))^2 +z1^2))));
    B2(i)=B2(i)+ ((1/pi)*r_drive*((2-k2square)*K2 - 2*E2)/(k2square*(sqrt((r_drive+r_ring(i))^2 +z2^2))));
    end
end
B=B1-B2;
% %----------------------------------------------------------------------------------------------------%    
% %----------------------3a. a_ij (depends on rings in film and trial value of λ)
% %----------------------------------------------------------------------------------------------------%
for i=1:N
    for j=1:N
        if j==i
A(i,j)= (lambda^2)/(d_eff*R)+ (delta_r*(log(16*pi*r_ring(j)/delta_r)-2))/(2*pi*R) ;
        else
            clear ksquare K E
       ksquare=4*r_ring(i)*r_ring(j)/(r_ring(i)+r_ring(j))^2;     
    [K,E]=ellipke(ksquare);
A(i,j)= ((1-ksquare/2)*K-E)*(r_ring(i)+r_ring(j))/(2*pi*N*r_ring(i));
        end
    end
end
% %----------------------------------------------------------------------------------------------------%    
% %----------------------3b. Solving a_ij c_j = b_i to get J_i
% %----------------------------------------------------------------------------------------------------%
C=A\B;
J=Id*(C./(d_eff*R));
%eval(sprintf('lambda%d =J',ls));
%----------------------------------------------------------------------------------------------------%    
%----------------------4. Vector potential at pickup coil (A due to film and drive coil)
%----------------------------------------------------------------------------------------------------%

Ap_d1=zeros(npr,nph);%Vector potential at j_th turn of i_th layer of pickup coil due to clockwise section of drive coil
Ap_d2=zeros(npr,nph);%Vector potential at j_th turn of i_th layer of pickup coil due to anti-clockwise section of drive coil
Ap_d=zeros(npr,nph);
for i=1:npr %layers in pickup coil
    for j=1:nph %turns in each layer of pickup coil
        for l=1:ndh %turns in each clockwise/anticlockwise section of drive coil
    clear k1square k2square z1 z2 K1 E1 K2 E2
        z1=sb+d+h_drive+h_pickup+((l-1)*wire)+((j-1)*wire);%distance of clockwise loop from pickup loop
        z2=sb+d+h_drive+h_pickup+((l+17)*wire)+((j-1)*wire);%distance of anti-clockwise loop from pickup loop
    k1square=(4*r_drive*(r_pickup+(i-1)*wire)/((r_drive+(r_pickup+(i-1)*wire))^2 +z1^2));
    k2square=(4*r_drive*(r_pickup+(i-1)*wire))/((r_drive+(r_pickup+(i-1)*wire))^2 +z2^2);
    [K1,E1]=ellipke(k1square);
    [K2,E2]=ellipke(k2square);
    Ap_d1(i,j)=Ap_d1(i,j)+ ((mu_o*(1/(4*pi))*4*Id*r_drive*((2-k1square)*K1 - 2*E1))/(k1square*(sqrt((r_drive+(r_pickup+(i-1)*wire))^2 +z1^2))));
    Ap_d2(i,j)=Ap_d2(i,j)+ ((mu_o*(1/(4*pi))*4*Id*r_drive*((2-k2square)*K2 - 2*E2))/(k2square*(sqrt((r_drive+(r_pickup+(i-1)*wire))^2 +z2^2))));
        end
    end
end
Ap_d=Ap_d1-Ap_d2;
%..............................%

Ap_r=zeros(npr,nph);%Vector potential at j_th turn of i_th layer of pickup coil due to rings in sample film
for i=1:npr %layers in pickup coil
    for j=1:nph %turns in each layer of pickup coil
        for l=1:N %number of rings in sample
            clear k1square z1 K1 E1 
        z1=h_pickup+((j-1)*wire);%distance of ring from pickup loop 
    k1square=((4*r_ring(l)*(r_pickup+(i-1)*wire))/(((r_ring(l)+(r_pickup+(i-1)*wire)))^2 +z1^2));
    [K1,E1]=ellipke(k1square);
    Ap_r(i,j)=Ap_r(i,j)+ ((mu_o*(1/(4*pi))*4*delta_r*d_eff*J(l)*r_ring(l)*((2-k1square)*K1 - 2*E1)/(k1square*(sqrt((r_ring(l)+(r_pickup+(i-1)*wire))^2 +z1^2)))));
        end
    end
end

%----------------------------------------------------------------------------------------------------%    
%----------------------5. Mutual inductance of pickup coil M= 2πr Ā/ I(due to film and drive coil)Ap_d(i,j)--Ap_r(i,j)
%----------------------------------------------------------------------------------------------------%
M=zeros(npr,nph);
for i=1:npr
    for j=1:nph
   M(i,j)=((Ap_d(i,j)-Ap_r(i,j))*(2*pi*(r_pickup+(i-1)*wire))/Id);%
    end
end
M_total(ls,1)=lambda;
M_total(ls,2)=sum(sum(M));
ls=ls+1
end
 %end

