clear;
clc;
for y=1:3; 
lam0=550*1e-9;
E0=8.85e-12; %Permittivity of free space
Mu0=12.5664e-7; %Permeability of free space
c=1/sqrt(E0*Mu0); %Speed of light (m/s)

Z(1)=1;
Mu(1)=1;
Ei(1)=1;

 Z(3)=1;
   Mu(3)=1;
    Ei(3)=1;
    
   display(2);
   Mu(2)=input('Enter Permeability(µ) = ');
   Ei(2)=input('Enter Permittivity(€) = ');
   
  
   t=(0:0.0001:1);
   t=t*550*1e-9;  
   
   if Ei(2)<0
 E=abs(Ei(2));
 else
 E=Ei(2);
   end
    lame=lam0/sqrt(E);

  if Mu(2)<0
 M=abs(Mu(2));
 else
 M=Mu(2);
  end
  
  for k =1:length(t)
      Beta1=(2*pi)./lame;
      total=1;
      w(k)=(c*2*pi)./lame;
      
      if Ei(2)>0 & Mu(2)>0 ;
      Z(2)=sqrt(M*E);
      end
       if Ei(2)<0 & Mu(2)<0 ;
      Z(2)=-sqrt(M*E);
       end
       if Ei(2)<0 & Mu(2)>0 ;
      Z(2)=j*sqrt(M*E);
      end
    
      for i=2:3
     ro(i)=(Z(i)-Z(i-1))/(Z(i)+Z(i-1));
     tao(i)=ro(i)+1;
      end
     
      if Ei(2)>0 & Mu(2)>0 ;
 A=exp(j*Beta1*t(k));
 Ac=exp(-j*Beta1*t(k));
 total=total*((1/tao(2))*[A,ro(2)*Ac;ro(2)*A,Ac]);
      end
  if Ei(2)<0 & Mu(2)<0;
    A=exp(-j*Beta1*t(k));
    Ac=exp(j*Beta1*t(k)); 
    total=total*((1/tao(2))*[A,ro(2)*Ac;ro(2)*A,Ac]);
  end
if Ei(2)<0 & Mu(2)>0;
     A=exp(Beta1*t(k));
    Ac=exp(-Beta1*t(k));
     total=total*((1/tao(2))*[A,ro(2)*Ac;ro(2)*A,Ac]);
      end

  
  total=1/tao(3)*total*[1 ro(3); ro(3) 1];
 rot(k)=total(2,1)/total(1,1);
%  taot(k)=1/total(1,1);
taot(k)=1-abs(rot(k));
  end
  
  if y==1;
      rot_1=rot;
      taot_1=taot;
      lame_1=lame;
  end
 
if y==2;
    taot_2=taot;
    rot_2=rot;
     lame_2=lame;
end
if y==3;
    taot_3=taot;
    rot_3=rot;
     lame_3=lame;
end
if y==4;
    taot_4=taot;
    rot_4=rot;
     lame_4=lame;
end
end

subplot(2,1,1);
plot(t/abs(lame_1),abs(rot_1),'.r');
hold on
plot(t/abs(lame_2),abs(rot_2),'.b');
hold on
plot(t/abs(lame_3),abs(rot_3),'.g');
legend('(E=-9]', '(E=-16]' ,'(E=-25]');
axis([0.01 1 0 1]);
grid;
xlabel('h/lambdae');
ylabel('Reflection');
title('Single DN Metamaterial Slab Length Vs Reflection');

subplot(2,1,2);
plot(t/abs(lame_1),abs(taot_1),'.r');
hold on
plot(t/abs(lame_2),abs(taot_2),'.b');
hold on
plot(t/abs(lame_3),abs(taot_3),'.g');
legend('(E=-9]', '(E=-16]' ,'(E=-25]');
axis([0.01 1 0 1]);
grid;
xlabel('h/lambdae');
ylabel('Transmission');
title('Single DN Metamaterial Slab Length Vs Transmission');