clear
%NEWTON RAPHSON METHOD
clc;
 
%Specified bus quantities
%Column 1: Bus, Column 2: Real Load demand, Column 3: Reactive load demand,
%Column 4: Real power generation, Column 5: Reactive power generation 
gbus = [1 2.0 1.0 0.0 0.0
2 0.0 0.0 0.5 1.0
3 1.5 0.6 0.0 0.0];
 
%Line Impedance
%Column 1: From, Column 2: To, Column 3: R pu, Column 4: X pu  
zdata=[1 2 0.02 0.08; 1 3 0.02 0.08; 2 3 0.02 0.08]
 
%ybus generation code from line impedance.
nl=zdata(:,1); nr=zdata(:,2); R=zdata(:,3); X=zdata(:,4); nbr=length(zdata(:,1)); nbus=max(max(nl),max(nr));
Z=R+j*X;
y=ones(nbr,1)./Z;
ybus=zeros(nbus,nbus);
 
% calculating off diagonal elements of YBus
for k=1:nbr;
    if nl(k)>0 & nr(k)>0
       ybus(nl(k),nr(k))=ybus(nl(k),nr(k))-y(k);
       ybus(nr(k),nl(k))=ybus(nl(k),nr(k));
    end
end
% calculating diagonal elements of YBus
for n=1:nbus
    for k=1:nbr
        if nl(k)==n|nr(k)==n
        ybus(n,n)=ybus(n,n)+y(k);
        else, end
    end
end
 
%Initializing values of voltages and tolerances:
t= 0.01;
v1=1.04+j*0;
v2=1+j*0;
v3=1.04+j*0;
 
del1=angle(v1);
del2=angle(v2);
del3=angle(v3);
 
%Loop for iterative calculation: 
for m=1:20
 
%Calculating P2, P3 (Real) and Q3 (Reactive power)
p2=(abs(v2)*abs(v1)*abs(ybus(2,1))*cos((angle(ybus(2,1)))+del1-del2))+abs(v2)*abs(v2)*abs(ybus(2,2))*cos(angle(ybus(2,2)))+(abs(v2)*abs(v3)*abs(ybus(2,3))*cos((angle(ybus(2,3)))+del3-del2));
q2=-(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))-abs(v2)*abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2))))-(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
p3=(abs(v3)*abs(v1)*abs(ybus(3,1))*cos((angle(ybus(3,1)))+del1-del3))+abs(v3)*abs(v3)*abs(ybus(3,3))*cos((angle(ybus(3,3))))+(abs(v2)*abs(v3)*abs(ybus(3,2))*cos((angle(ybus(3,2)))+del2-del3));
 
%Changes in P and Q = P(Specified)-P(Calculated). 
%P(Specified) = Power Generation-Load Demand.
delp20=gbus(2,4)-gbus(2,2)-p2;
delp30=gbus(3,4)-gbus(3,2)-p3;
delq20=gbus(2,5)-gbus(2,3)-q2;
 
%Reactive power at bus 3 calculated using
q3=-(abs(v3)*abs(v1)*abs(ybus(3,1))*sin((angle(ybus(3,1)))+del1-del3))-abs(v3)*abs(v3)*abs(ybus(3,3))*sin((angle(ybus(3,3))))-(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
 
%Reactive power generation at bus 3:
qg3=gbus(3,3)+q3;
 
%Checking constrainet on Reactive power generation at bus 3:
if qg3>=0 & qg3<=1.5
    
%Jacobian Matrix Generation 
J(1,1)=(abs(v2)*abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))+(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
J(1,2)=-(abs(v2)*abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
J(1,3)=(abs(v1)*abs(ybus(2,1))*cos((angle(ybus(2,1)))+del1-del2))+2*(abs(v2)*abs(ybus(2,2))*cos((angle(ybus(2,2)))))+(abs(v3)*abs(ybus(2,3))*cos((angle(ybus(2,3)))+del3-del2));
J(2,1)=-(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
J(2,2)=(abs(v3)*abs(v1)*abs(ybus(3,1))*sin((angle(ybus(3,1)))+del1-del3))+(abs(v3)*abs(v2)*abs(ybus(3,2))*sin((angle(ybus(3,2)))+del2-del3));
J(2,3)=(abs(v3)*abs(ybus(3,2))*cos((angle(ybus(3,2)))+del2-del3));
J(3,1)=(abs(v2)*abs(v1)*abs(ybus(2,1))*cos((angle(ybus(2,1)))+del1-del2))+(abs(v2)*abs(v3)*abs(ybus(2,3))*cos((angle(ybus(2,3)))+del3-del2));
J(3,2)=-(abs(v2)*abs(v3)*abs(ybus(2,3))*cos((angle(ybus(2,3)))+del3-del2));
J(3,3)=-(abs(v1)*abs(ybus(2,1))*sin((angle(ybus(2,1)))+del1-del2))-2*(abs(v2)*abs(ybus(2,2))*sin((angle(ybus(2,2)))))-(abs(v3)*abs(ybus(2,3))*sin((angle(ybus(2,3)))+del3-del2));
 
%Calculating bus 2 & 3 voltage angle and bus 2 voltage magnitude:
inv(J);
A=[del2;del3;abs(v2)];
delA0=[delp20;delp30;delq20];
delA1=inv(J)*delA0;
delA1;
b0=abs(v2);
A1=[del2;del3;b0]+delA1;
A1-delA0;
 
%Cheking power missmatch tolerance 
    if(A1-A)<=t & max(delA0)<=0.01
        break;
end
A1
%Updated value after each iteration.
del2=A1(1,1);
del3=A1(2,1);
v2=A1(3,1);
 
else
    disp('PV bus')
end
 
end
 
%Aparent power calculation in all buses. 
p1=(abs(v1)*abs(v2)*abs(ybus(1,2))*cos((angle(ybus(1,2)))+del2-del1))+abs(v1)*abs(v1)*abs(ybus(1,1))*cos(angle(ybus(1,1)))+(abs(v1)*abs(v3)*abs(ybus(1,3))*cos((angle(ybus(1,3)))+del3-del1));
q1=-(abs(v1)*abs(v2)*abs(ybus(1,2))*sin((angle(ybus(1,2)))+del2-del1))-abs(v1)*abs(v1)*abs(ybus(1,1))*sin((angle(ybus(1,1))))-(abs(v1)*abs(v3)*abs(ybus(1,3))*sin((angle(ybus(1,3)))+del3-del1));
S1=p1+j*q1;
S2=p2+j*q2;
S3=p3+j*q3;
 
[d,e]=pol2cart(del2,v2);
v2c=d+j*e;
[g,h]=pol2cart(del3,v3);
v3c=g+j*h;
 
SL=zeros(3,3);
 
I12=-ybus(1,2)*(v1-v2c);
I21=-I12;
I13=-ybus(1,3)*(v1-v3c);
I31=-I13;
I23=-ybus(2,3)*(v2c-v3c);
I32=-I23;
 
SL(1,2)=v1*conj(I12);
SL(2,1)=v2c*conj(I21);
SL(1,3)=v1*conj(I13);
SL(3,1)=v3c*conj(I31);
SL(2,3)=v2c*conj(I23);
SL(3,2)=v3c*conj(I32);
 
disp('Ybus Matrix:')
ybus
disp('Solution converges in iteration:')
m
disp('Reactive power generation at bus 3:'); 
qg3
disp('Angle of voltage at bus 2:')
del2
disp('Angle of voltage at bus 3:')
del3
disp('Bus 2 voltage magnitude:')
v2
disp('Aparent power at bus 1,2 &3')
 
S1=p1+j*q1
S2=p2+j*q2
S3=p3+j*q3
 
disp('Line flows')
SL
 
%Tranmission loss calculation:
Loss1=SL(1,2)+ SL(2,1);
Loss2=SL(1,3)+ SL(3,1);
Loss3=SL(2,3)+ SL(3,2);
 
TransmissionLoss=real(Loss1)+real(Loss2)+real(Loss3)
