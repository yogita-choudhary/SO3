N=15001;
eo=zeros(3,N+1);
ew=zeros(3,N+1);
mo=zeros(3,N+1);
mw=zeros(3,N+1);
euw=zeros(3,N+1);
s=zeros(3,N+1);
Tc=zeros(3,N+1);
H=zeros(3,N+1);
R=zeros(3,3,N+1);
e1=[1;0;0];
e2=[0;1;0];
e3=[0;0;1];
for agf=1:N+1
   R(:,:,agf)=eye(3);
end
R(:,:,1)=expm(Hat(e3)*(-pi/4))*expm(Hat(e2)*pi/6)*expm(Hat(e1)*(pi/3));
temp=zeros(1,N+1);
integral=[0;0;0];
omega=zeros(3,N+1);
omega(:,1)=[-0.12;0.42;-0.38];
E=zeros(3,3,N+1);
eta1=2.6;
eta2=2.6;
eta3=2.6;
eta4=2.6;
eta5=2.6;
eta6=2.6;
tf=5;
B=zeros(1,N+1);
dt=0.001;
syms t;
v=zeros(3,N+1);
k1=5.15;
k2=1;
k3=15;
k4=14;
v(:,1)=[0;0;0];
[omegad,Rd]=desiredw_tr();
g=[0;0;0];
k=zeros(3,N+1);
J=[3,0,0;0,2,0;0,0,1];
f=0;
unom=[0;0;0];
udisc=[0;0;0];
eo(:,1)=(1/(2*sqrt(1+trace(Rd(:,:,1)'*R(:,:,1)))))*(Vee(Rd(:,:,1)'*R(:,:,1)-R(:,:,1)'*Rd(:,:,1)));
euw(:,1)=omega(:,1)-R(:,:,1)'*Rd(:,:,1)*omegad(:,1);
RdtR(:,:,1)=Rd(:,:,1)'*R(:,:,1);
E(:,:,1)=(1/(2*sqrt(1+trace(RdtR(:,:,1)))))*(trace(RdtR(:,:,1)')*eye(3)-RdtR(:,:,1)'+2*eo(:,1)*eo(:,1)');
ew(:,1)=E(:,:,1)*euw(:,1);
etemp=zeros(1,N+1);
B(1,1)=2-sqrt(1+trace(Rd(:,:,1)'*R(:,:,1)));
for n=1:N
    d(:,n)=[2*sin(4*pi*n*dt);2*sin(4*pi*n*dt);2*sin(4*pi*n*dt)];
    etemp(1,n)=trace(Rd(:,:,n)'*R(:,:,n));
    RdtR(:,:,n)=Rd(:,:,n)'*R(:,:,n);
    euw(:,n)=omega(:,n)-R(:,:,n)'*Rd(:,:,n)*omegad(:,n);
    E(:,:,n)=(1/(2*sqrt(1+trace(RdtR(:,:,n)))))*(trace(RdtR(:,:,n)')*eye(3)-RdtR(:,:,n)'+2*eo(:,n)*eo(:,n)');
    invE=inv(E(:,:,n));
    R(:,:,n+1)=R(:,:,n)*expm(Hat(omega(:,n))*dt);
    tg=tf-n*dt;
    RdtR(:,:,n+1)=Rd(:,:,n+1)'*R(:,:,n+1);
    phi(1)=eta1*(exp(eo(1,n))-1)/(exp(eo(1,n))*(tg));
    phi(2)=eta2*(exp(eo(2,n))-1)/(exp(eo(2,n))*(tg));
    phi(3)=eta3*(exp(eo(3,n))-1)/(exp(eo(3,n))*(tg));
    phidot(1)=eta1*(exp(eo(1,n))-1)/(exp(eo(1,n))*(tg)^2);
    phidot(2)=eta2*(exp(eo(2,n))-1)/(exp(eo(2,n))*(tg)^2);
    phidot(3)=eta3*(exp(eo(3,n))-1)/(exp(eo(3,n))*(tg)^2);
    phieo(1)=eta1*(exp(-eo(1,n)))/(tg);
    phieo(2)=eta2*(exp(-eo(2,n)))/(tg);
    phieo(3)=eta3*(exp(-eo(3,n)))/(tg);
    beta(1)=eta4*(exp(ew(1,n)+phi(1))-1)/(exp(ew(1,n)+phi(1))*(tg));
    beta(2)=eta5*(exp(ew(2,n)+phi(2))-1)/(exp(ew(2,n)+phi(2))*(tg));
    beta(3)=eta6*(exp(ew(3,n)+phi(3))-1)/(exp(ew(3,n)+phi(3))*(tg));
    eo(:,n+1)=eo(:,n)+ew(:,n)*dt;
    E(:,:,n+1)=(1/(2*sqrt(1+trace(RdtR(:,:,n+1)))))*(trace(RdtR(:,:,n+1)')*eye(3)-RdtR(:,:,n+1)'+2*eo(:,n+1)*eo(:,n+1)');
    alphad=-Hat(omega(:,n))*R(:,:,n)'*Rd(:,:,n)*omegad(:,n)+R(:,:,n)'*Rd(:,:,n)*((omegad(:,n+1)-omegad(:,n))/dt);
    if n*dt>tf-0.005
       s(:,n)=ew(:,n)-ew(:,1)-integral;
      udisc=(-k1*s(:,n)/sqrt(sqrt(s(1,n)^2+s(2,n)^2+s(3,n)^2)))+v(:,n)-k2*s(:,n);
      v(:,n+1)=(-k3*(s(:,n)/sqrt(s(1,n)^2+s(2,n)^2+s(3,n)^2))-k4*s(:,n))*dt+v(:,n);
      unom=[0;0;0];
    else
        unom=-beta'-eo(:,n)-phidot'-[phieo(1) 0 0;0 phieo(2) 0; 0 0 phieo(3)]*ew(:,n);
        integral=integral+unom*dt;
        s(:,n)=ew(:,n)-ew(:,1)-integral;
        udisc=(-k1*s(:,n)/sqrt(sqrt(s(1,n)^2+s(2,n)^2+s(3,n)^2)))+v(:,n)-k2*s(:,n);
        v(:,n+1)=(-k3*(s(:,n)/sqrt(s(1,n)^2+s(2,n)^2+s(3,n)^2))-k4*s(:,n))*dt+v(:,n);
    end
    g=unom+udisc;
    wjw(:,n)=Hat(omega(:,n))*J*omega(:,n);
    Tc(:,n)=J*invE*(g+E(:,:,n)*alphad-(E(:,:,n+1)-E(:,:,n))*(omega(:,n)-R(:,:,n)'*Rd(:,:,n)*omegad(:,n))/dt)+wjw(:,n);
    ew(:,n+1)=ew(:,n)+(g+E(:,:,n)*inv(J)*d(:,n))*dt;
    omega(:,n+1)=omega(:,n)+(inv(J)*(-wjw(:,n)+Tc(:,n)+d(:,n)))*dt;
    B(n)=real(2-sqrt(1+trace(Rd(:,:,n)'*R(:,:,n))));
    H(:,n)=E(:,:,n)*inv(J)*d(:,n);
end
gem=0.00:0.001:15.0010;
figure(1)
plot(gem,omega(1,1:N+1),'DisplayName','actual omega1')
hold on 
plot(gem,omegad(1,1:N+1),'DisplayName','Desired omega1')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
legend('{\Omega_{1}}','{\Omega_{desired}}')
xlabel('Time (s)','interpreter','latex') 
ylabel('Angular velocity (rad/s)','interpreter','latex') 

figure(2)
plot(gem,omega(2,1:N+1),'DisplayName','actual omega2')
hold on 
plot(gem,omegad(2,1:N+1),'DisplayName','Desired omega2')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
legend('{\Omega_{2}}','{\Omega_{desired}}')
xlabel('Time (s)','interpreter','latex') 
ylabel('Angular velocity (rad/s)','interpreter','latex') 

figure(3)
plot(gem,omega(3,1:N+1),'DisplayName','actual omega3')
hold on 
plot(gem,omegad(3,1:N+1),'DisplayName','Desired omega3')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
legend('{\Omega_{3}}','{\Omega_{desired}}')
xlabel('Time (s)','interpreter','latex') 
ylabel('Angular velocity (rad/s)','interpreter','latex') 

figure(4)
plot(gem,B,'DisplayName','psii')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
ylabel('Attitude error function (${\varphi}$)','interpreter','latex')
xlabel('Time (s)')


figure(5)
plot(gem,Tc(1,1:N+1),'DisplayName','\tau_{c_1}')
hold on 
plot(gem,Tc(2,1:N+1),'DisplayName','\tau_{c_2}')
hold on 
plot(gem,Tc(3,1:N+1),'DisplayName','\tau_{c_3}')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
legend('{u_1}','{u_2}','{u_3}','interpreter','latex')
xlabel('Time (s)','interpreter','latex') 
ylabel('Control torque (Nm)','interpreter','latex') 

figure(6)
plot(gem,euw(1,1:N+1),'DisplayName','\tau_{c_1}')
hold on 
plot(gem,euw(2,1:N+1),'DisplayName','\tau_{c_2}')
hold on 
plot(gem,euw(3,1:N+1),'DisplayName','\tau_{c_3}')
hold on
plot([1 1]*5, ylim, '--k') 
hold off
legend('e_{v_1}','e_{v_2}','e_{v_3}','interpreter','latex')
xlabel('Time (s)','interpreter','latex') 
ylabel('Angular velocity error (rad/s)','interpreter','latex') 
toc;