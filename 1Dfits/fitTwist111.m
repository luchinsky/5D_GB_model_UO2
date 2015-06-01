function params=fitTwist111(twist111,params)
ergb=params(1);
ang=pi/3*180/pi;
C=twist111(abs((twist111(:,1)-ang))==min(abs(twist111(:,1)-ang)),2);
a=[params(38) params(37) params(36)];
a(1)=a(1)*ergb;

da=0.01;
theta=twist111(:,1)';
E=twist111(:,2)';
[theta, i]=sort(theta);
E=E(i);

A=stepParameters(C,a,da,theta,E);
A=[C A];

dx=(max(theta)-min(theta))/200;
x=min(theta):dx:max(theta);
y=equation(x,A);
h=figure(1);
plot(theta,E,'.',x,y,'-')

A(1)=A(1)/ergb;
A(2)=A(2)/ergb;

title(strcat('UO2 111 Twist   eRGB: ',num2str(ergb),' Jm^{-2}'))
text(20,0.9,strcat('Max Energy: ', num2str(A(2))))
text(20,0.75,strcat('\epsilon_{twin}: ', num2str(A(1))))
text(20,0.6,strcat('Max Angle: ', num2str(A(3))))
text(20,0.45,strcat('Shaping Factor: ', num2str(A(4))))
xlabel('Angle')
ylabel('Energy')
saveas(h,'./Results/111Twist.jpg')

params(36)=A(4);
params(37)=A(3);
params(38)=A(2);
params(28)=A(1);

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
        title('111 Twist')
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
theta3=pi/3;

en2=a(2);
en3=a(1);
theta2=a(3);
a1=a(4);

a2=a1;

cx1=x>=theta1;
cx2=x<theta2;
cx=cx1&cx2;
dtheta=theta2-theta1;
theta=x(cx);
thetan=(theta-theta1)./dtheta*pi/2;
E1=en1+(en2-en1).*(sin(thetan)-a1.*(sin(thetan).*log(sin(thetan))));
E1(thetan==0)=en1;

cx1=x>=theta2;
cx2=x<=theta3;
cx=cx1&cx2;
dtheta=theta2-theta3;
theta=x(cx);
thetan=(theta-theta3)./dtheta*pi/2;
E2=en3+(en2-en3).*(sin(thetan)-a2.*(sin(thetan).*log(sin(thetan))));
E2(thetan==0)=en3;

E=[E1 E2];
end

function X2=leastSquare(g,y)

X2=sum((y-g).^2);
end