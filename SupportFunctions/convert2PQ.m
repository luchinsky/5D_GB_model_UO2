function [P,Q]=convert2PQ(data)
    R=[cos(pi/2) 0 sin(pi/2); 0 1 0; -sin(pi/2) 0 cos(pi/2)];
   
    x1=[data(1) data(2) data(3)];
    y1=[data(4) data(5) data(6)];
    z1=[data(7) data(8) data(9)];
    x2=[data(10) data(11) data(12)];
    y2=[data(13) data(14) data(15)];
    z2=[data(16) data(17) data(18)];
    
    P(1,:)=x1/norm(x1);
    P(2,:)=y1/norm(y1);
    P(3,:)=z1/norm(z1);
    
    Q(1,:)=x2/norm(x2);
    Q(2,:)=y2/norm(y2);
    Q(3,:)=z2/norm(z2);
    
    P=R*P;
    Q=R*Q;
end