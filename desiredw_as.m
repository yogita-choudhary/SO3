function [omegad,Rdd] = desiredw_as()
syms t;
fi=0;
si=0;
ki=0;
e1=[1;0;0];
e2=[0;1;0];
e3=[0;0;1];
Rd=expm(Hat(e3)*ki)*expm(Hat(e2)*si)*expm(Hat(e1)*fi);
Rd_dot=diff(Rd,t);
omhat=Rd'*Rd_dot;
N=15;
omegad=zeros(3,15002);
omegahat=zeros(3,3);
for agf=1:15002
    Rdd(:,:,agf)=eye(3);
end
k=1;
for p=0.00:0.001:N+0.0010
    omegahat=real(subs(omhat,t,p));
    Rdd(:,:,k)=real(subs(Rd,t,p));
    omegad(1,k)=real(omegahat(3,2));
    omegad(2,k)=real(omegahat(1,3));
    omegad(3,k)=real(omegahat(2,1));
    k=k+1;
    
end
end