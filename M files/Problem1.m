close all;
clear all;
clc;
zdata=[1 2 0.05 0.15; 1 3 0.1 0.3; 2 3 0.15 0.45; 2 4 0.1 0.3; 3 4 0.05 0.15]; %Data Input
nl=zdata(:,1); nr=zdata(:,2); R=zdata(:,3); X=zdata(:,4); nbr=length(zdata(:,1)); nbus=max(max(nl),max(nr));
Z=R+X*i;
y=ones(nbr,1)./Z;
Y=zeros(nbus,nbus);
% calculating off diagonal elements of YBus
for k=1:nbr
    if nl(k)>0 && nr(k)>0
       Y(nl(k),nr(k))=Y(nl(k),nr(k))-y(k);
       Y(nr(k),nl(k))=Y(nl(k),nr(k));
    end
end
% calculating diagonal elements of YBus
for n=1:nbus
    for k=1:nbr
    if nl(k)==n|nr(k)==n
    Y(n,n)=Y(n,n)+y(k);
    else, end
    end
end
YBUS=Y %YBUS Matix
%Given Real and Reactive power of buses:
p2=0.5;
q2=-0.2;
p3=-1;
q3=0.5;
p4=0.3;
q4=-0.1;
 
v1=1.04; %Slack bus voltage
%initializing voltages
v2=1-0i;
v3=1-0i;
v4=1-0i;
dx1=1+0i;
dx2=1+0i;
dx3=1+0i;
% if bus 2 is PQ bus Calculation of voltages:
o=1;
for k=1:o
    g2=(((p2-i*q2)/conj(v2))+(-Y(1,2)*v1)+(-Y(2,3)*v3)+(-Y(2,4)*v4))/Y(2,2);
    g3=(((p3-i*q3)/conj(v3))+(-Y(1,3)*v1)+(-Y(2,3)*g2)+(-Y(3,4)*v4))/Y(3,3);
    g4=(((p4-i*q4)/conj(v4))+(-Y(1,4)*v1)+(-Y(2,4)*g2)+(-Y(3,4)*g3))/Y(4,4);
    dx1=g2-v2;
    dx2=g3-v3;
    dx3=g4-v4;
    if abs(dx1)&&abs(dx2)&&abs(dx3)<=0.00001
        break
    end
    v2=v2+dx1;
    v3=v3+dx2;
    v4=v4+dx3;
end
value1=[v2;v3;v4]
 
%let bus 2 is PV bus instead of PQ
v2=1.04+0i;   
%initializing voltage
v3=1-0i;
v4=1-0i;
dx1=1+0i;
dx2=1+0i;
dx3=1+0i;
 
%Calculation of voltages and Reactive power:
for k=1:o    % if 0.20<Q<1 tolerance  Condition for checking PV bus
    q2=-imag(conj(v2)*((Y(2,1)*v1)+(Y(2,2)*v2)+(Y(2,3)*v3)+(Y(2,4)*v4))); % calculating the value of PV bus Q2  
    if q2<=0.20 && q2>1
        break
    end
    del=((p2-q2*i)/conj(v2)-Y(1,2)*v1-Y(2,3)*v3-Y(2,4)*v4)/Y(2,2);
    del2=angle(del); %for matlab this angle is in radian
    g2=1.04*(cos(del2)+sin(del2)*i);
    g3=(((p3-i*q3)/conj(v3))+(-Y(1,3)*v1)+(-Y(2,3)*g2)+(-Y(3,4)*v4))/Y(3,3);
    g4=(((p4-i*q4)/conj(v4))+(-Y(1,4)*v1)+(-Y(2,4)*g2)+(-Y(3,4)*g3))/Y(4,4);
    dx1=g2-v2;
    dx2=g3-v3;
    dx3=g4-v4;
    if abs(dx1)&&abs(dx2)&&abs(dx3)<=0.00001
        break
    end
    v2=v2+dx1;
    v3=v3+dx2;
    v4=v4+dx3;
end
value2=[v2;v3;v4]

    v2=1.04+j*0;
    v3=1-0i;
    v4=1-0i;
    dx1=1+0i;
    dx2=1+0i;
    dx3=1+0i; 
 %Calculation of voltages and Reactive power: 
for k=1:o  % if 0.25<Q<1 tolerance
    q2=-imag(conj(v2)*((Y(2,1)*v1)+(Y(2,2)*v2)+(Y(2,3)*v3)+(Y(2,4)*v4))); % calculating the value of PV bus Q2
    if q2>=0.25 & q2<=1.0 % Condition for checking PV bus
    del=((p2-q2*i)/conj(v2)-Y(1,2)*v1-Y(2,3)*v3-Y(2,4)*v4)/Y(2,2);
    del2=angle(del); %for matlab this angle is in radian
    g2=1.04*(cos(del2)+sin(del2)*i);
    %Outside this range bus 2 becomes PQ bus and q2=q2(minimum)and v2=1+j*0.
    else
    q2=0.25;
    v2=1+j*0;
    g2=(((p2-i*q2)/conj(v2))+(-Y(1,2)*v1)+(-Y(2,3)*v3)+(-Y(2,4)*v4))/Y(2,2);
    end
      g3=(((p3-i*q3)/conj(v3))+(-Y(1,3)*v1)+(-Y(2,3)*g2)+(-Y(3,4)*v4))/Y(3,3);
    g4=(((p4-i*q4)/conj(v4))+(-Y(1,4)*v1)+(-Y(2,4)*g2)+(-Y(3,4)*g3))/Y(4,4);
    dx1=g2-v2;
    dx2=g3-v3;
    dx3=g4-v4;
    if abs(dx1)&&abs(dx2)&&abs(dx3)<=0.00001
        break
    end
    v2=v2+dx1;
    v3=v3+dx2;
    v4=v4+dx3;
end
value3=[v2;v3;v4]
