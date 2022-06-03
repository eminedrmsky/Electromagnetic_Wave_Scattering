clear;
clc;
for y=1:3;
N = input('Enter the number of layers= ');

lam0=(400:0.01:700)*1e-9;

E0=8.85e-12; %Permittivity of free space
Mu0=12.5664e-7; %Permeability of free space
c=1/sqrt(E0*Mu0); %Speed of light (m/s)

for layer=1:N ;
    if layer==1;
        Z(layer)=1;
        Mu(layer)=1;
        Ei(layer)=1;
        t(layer)=100000000000000000000000000000000000;
        
        elseif layer==N;
            Z(layer)=1.50;
           Mu(layer)=1;
           Ei(layer)=1;
           t(layer)=100000000000000000000000000000000000;
           
        else   
           display(layer);
           Mu(layer)=input('Enter Permeability(µ) = ');
           Ei(layer)=input('Enter Permittivity(€) = ');
           t(layer)=input('Enter Layer Length (m) = ');
       %    Z(layer)=input('enter index');
       end
end
t=t*550*1e-9;    
for i=1:N
 if Ei(i)<0
 E(i)=abs(Ei(i));
 else
 E(i)=Ei(i);
 end
end

for i=1:N
 if Mu(i)<0
 M(i)=abs(Mu(i));
 else
 M(i)=Mu(i);
 end
end

for k =1:length(lam0) 

 Beta1=(2*pi)./lam0(k);
 total=1;
 w(k)=(c*2*pi)./lam0(k);
 
 for i=1:N
     if Ei(i)>0 & Mu(i)>0;
       Z(i)=sqrt(M(i)*E(i));  
     end
     if Ei(i)<0 & Mu(i)>0;
       Z(i)=j*sqrt(E(i)*M(i));
     end
     if Ei(i)<0 & Mu(i)<0;
         Z(i)=-sqrt(E(i)*M(i));
     end
 end
 
  ro(1)=(Z(1)-1)/(Z(1)+1);
 tao(1)=ro(1)+1;
 
 for i=2:N
 ro(i)=(Z(i)-Z(i-1))/(Z(i)+Z(i-1));
 tao(i)=ro(i)+1;
 end
  for i =2:(N-1)
%  if Ei(i)<0
%  sign=-1;
%  else
%  sign=1;
 
%   end
%  Beta=Beta1*sign;
if Ei(i)>0 & Mu(i)>0 ;
 A=exp(j*Beta1*t(i)*sqrt(E(i)*M(i)));
 Ac=exp(-j*Beta1*t(i)*sqrt(E(i)*M(i)));
 total=total*((1/tao(i))*[A,ro(i)*Ac;ro(i)*A,Ac]);
 
elseif Ei(i)<0 & Mu(i)<0;
    A=exp(-j*Beta1*t(i)*sqrt(E(i)*M(i)));
    Ac=exp(j*Beta1*t(i)*sqrt(E(i)*M(i))); 
    total=total*((1/tao(i))*[A,ro(i)*Ac;ro(i)*A,Ac]);
    
else 
     A=exp(Beta1*t(i)*sqrt(E(i)*M(i)));
    Ac=exp(-Beta1*t(i)*sqrt(E(i)*M(i)));
     total=total*((1/tao(i))*[A,ro(i)*Ac;ro(i)*A,Ac]);
end
  end
  
  total=1/tao(N)*total*[1 ro(N); ro(N) 1];
 rot(k)=total(2,1)/total(1,1);
%  taot(k)=1/total(1,1);
taot(k)=1-abs(rot(k));
end
if y==1;
    taot_1=taot;
    rot_1=rot;
end
if y==2;
    taot_2=taot;
    rot_2=rot;
end
if y==3;
    taot_3=taot;
    rot_3=rot;
end
end
subplot(2,1,1);
plot(lam0,abs(rot_1),'.r');
hold on
plot(lam0,abs(rot_2),'.b');
hold on
plot(lam0,abs(rot_3),'.g');
grid;
xlabel('wavelength');
ylabel('Reflection');
legend('(E=-9]','(E=-16]', 'E=-25]' );
title('Half Wave Length Single DN Metamaterial Slab');

subplot(2,1,2);
plot(lam0,abs(taot_1),'.r');
hold on
plot(lam0,abs(taot_2),'.b');
hold on
plot(lam0,abs(taot_3),'.g');
grid;
xlabel('wavelength');
ylabel('Transmission');
legend('(E=-9]','(E=-16]', '(E=-25]');
title('Half Wave Length Single DN Metamatarial Slab');