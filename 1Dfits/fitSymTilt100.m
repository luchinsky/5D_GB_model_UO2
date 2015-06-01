function params=fitSymTilt100(stgb100,params)
ergb=params(1);
ang=acos(4/5)*180/pi;
c(1)=stgb100(abs((stgb100(:,1)-ang))==min(abs(stgb100(:,1)-ang)),2);
ang=acos(3/5)*180/pi;
c(2)=stgb100(abs((stgb100(:,1)-ang))==min(abs(stgb100(:,1)-ang)),2);
C=[c(1) c(2)];
a=[params(17) params(13) params(15) params(18) params(19)];
a(1)=a(1)*ergb;
a(2)=a(2)*ergb;
a(3)=a(3)*ergb;

da=0.1;
theta=stgb100(:,1)';
E=stgb100(:,2)';
[theta, i]=sort(theta);
E=E(i);

A=stepParameters(C,a,da,theta,E);
A=[C A];

dx=(max(theta)-min(theta))/200;
x=min(theta):dx:max(theta);
y=equation(x,A);

h=figure(1);
plot(theta,E,'.',x,y,'-');

A(1)=A(1)/ergb;
A(2)=A(2)/ergb;
A(3)=A(3)/ergb;
A(4)=A(4)/ergb;
A(5)=A(5)/ergb;

title(strcat('UO2 100 Symmetric Tilt   eRGB: ',num2str(ergb),' Jm^{-2}'))
text(20,1.2,strcat('\Sigma5_1: ', num2str(A(1))))
text(20,1.1,strcat('\Sigma5_2: ', num2str(A(2))))
text(20,1.0,strcat('\Sigma17: ', num2str(A(3))))
text(20,0.9,strcat('Max 1 Energy: ', num2str(A(4))))
text(20,0.8,strcat('Max 2 Energy: ', num2str(A(5))))
text(20,0.7,strcat('Max 1 Angle: ', num2str(A(6))))
text(20,0.6,strcat('Max 2 Angle: ', num2str(A(7))))
xlabel('Angle')
ylabel('Energy')
saveas(h,'./Results/100SymTilt.jpg')

params(13)=A(4);
params(14)=A(1);
params(15)=A(5);
params(16)=A(2);
params(17)=A(3);
params(18)=A(6);
params(19)=A(7);

end

function a=stepParameters(C,a,da,x,y)
tol=0.000000001;
it=0;
while da>tol
    it=it+1;
    g=equation(x,[C a]);
    now=leastSquare(g,y);
    for i=1:length(a)
        a(i)=a(i)+da;
        g=equation(x,[C a]);
        X2(i,1)=leastSquare(g,y);
        a(i)=a(i)-da;
    end
    for i=1:length(a)
        a(i)=a(i)-da;
        g=equation(x,[C a]);
        X2(i,2)=leastSquare(g,y);
        a(i)=a(i)+da;
    end
    [val ind]=min(X2(:));
    col=1;
    if ind>length(a)
        ind=ind-length(a);
        col=2;
    end
    if now<=min(X2(:))
        da=da/2;
    elseif col==1
        a(ind)=a(ind)+da;
    elseif col==2
        a(ind)=a(ind)-da;
    end
    
    if mod(it,5)==0
        dx=(max(x)-min(x))/200;
        px=min(x):dx:max(x);
        py=equation(px,[C a]);
        figure(1)
        plot(x,y,'.',px,py,'-')
        title('100 Symmetric Tilt')
        xlabel('Misorientation Angle')
        ylabel('Energy')
        drawnow
    end
end


end

function E=equation(x,a)
x=x*pi/180;

en1=0;
theta1=0;
theta3=acos(4/5);
theta5=acos(3/5);
theta7=pi/2;
theta6=2*acos(5/sqrt(34));
en7=en1;

en3=a(1);
en5=a(2);
en6=a(3);
en2=a(4);
en4=a(5);
theta2=a(6);
theta4=a(7);

a1=0.5;
a2=0.5;
a3=0.5;
a4=0.5;
a5=0.5;
a6=0.5;

cx1=x>=theta1;
cx2=x<theta2;
cx=cx1&cx2;
dtheta=theta2-theta1;
theta=x(cx);
thetan=(theta-theta1)./dtheta*pi/2;
E1=en1+(en2-en1).*(sin(thetan)-a1.*(sin(thetan).*log(sin(thetan))));
E1(thetan==0)=en1;

cx1=x>=theta2;
cx2=x<theta3;
cx=cx1&cx2;
dtheta=theta2-theta3;
theta=x(cx);
thetan=(theta-theta3)./dtheta*pi/2;
E2=en3+(en2-en3).*(sin(thetan)-a2.*(sin(thetan).*log(sin(thetan))));
E2(thetan==0)=en3;

cx1=x>=theta3;
cx2=x<theta4;
cx=cx1&cx2;
dtheta=theta4-theta3;
theta=x(cx);
thetan=(theta-theta3)./dtheta*pi/2;
E3=en3+(en4-en3).*(sin(thetan)-a3.*(sin(thetan).*log(sin(thetan))));
E3(thetan==0)=en3;

cx1=x>=theta4;
cx2=x<theta5;
cx=cx1&cx2;
dtheta=theta4-theta5;
theta=x(cx);
thetan=(theta-theta5)./dtheta*pi/2;
E4=en5+(en4-en5).*(sin(thetan)-a4.*(sin(thetan).*log(sin(thetan))));
E4(thetan==0)=en5;

cx1=x>=theta5;
cx2=x<theta6;
cx=cx1&cx2;
dtheta=theta6-theta5;
theta=x(cx);
thetan=(theta-theta5)./dtheta*pi/2;
E5=en5+(en6-en5).*(sin(thetan)-a5.*(sin(thetan).*log(sin(thetan))));
E5(thetan==0)=en5;

cx1=x>=theta6;
cx2=x<theta7;
cx=cx1&cx2;
dtheta=theta6-theta7;
theta=x(cx);
thetan=(theta-theta7)./dtheta*pi/2;
E6=en7+(en6-en7).*(sin(thetan)-a6.*(sin(thetan).*log(sin(thetan))));
E6(thetan==0)=en7;

E=[E1 E2 E3 E4 E5 E6];
end

function X2=leastSquare(g,y)
X2=sum((y-g).^2);
end