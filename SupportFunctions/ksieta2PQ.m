
% axis=[1 1 1];   %rotation set [1 0 0], [1 1 0], or [1 1 1]
% ksi=4;          %ksi angle in degrees
% eta=39;         %eta angle in degrees
% phi=pi/2;       %phi angle in radians 0 for twist, pi/2 for tilt
% [P Q]=ksietaPQ(axis, ksi, eta, phi)

function [P Q]=ksieta2PQ(axis, ksi, eta, phi)
%convert to radians
ksi=ksi*pi/180;
eta=eta*pi/180;

%calculate theta1 and theta2
theta1=0.5*(eta-ksi);
theta2=0.5*(eta+ksi);

%calculate rotation matrix for initial orientation
if axis==[1 0 0]
    R=eye(3);
else
    R=rotateVector(axis);
end

%check for [1 1 1] rotations and adjust to line up reference axis
if axis==[1 1 1]
    dir1=[1 -1 0]/sqrt(2);
    dir2=cross(axis/sqrt(3),dir1);
    dir=R(1,:);
    dtheta=-atan2(dir2*dir',dir1*dir');
else
    dtheta=0;
end

%apply ksi and eta rotations
P=[1 0 0; 0 cos(theta1+dtheta) -sin(theta1+dtheta); 0 sin(theta1+dtheta) cos(theta1+dtheta)]*R;
Q=[1 0 0; 0 cos(theta2+dtheta) -sin(theta2+dtheta); 0 sin(theta2+dtheta) cos(theta2+dtheta)]*R;

%apply phi rotations
P=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)]*P;
Q=[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)]*Q;

end

function R=rotateVector(v1)
v2=[1 0 0]; %set destination vector

v1=v1/norm(v1); %normalize vector

x=cross(v1,v2); %get rotation axis
x=x/norm(x);

theta=acos(dot(v1,v2)); %get rotation angle
A=[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

R=eye(3)+sin(theta)*A+(1-cos(theta))*A^2;   %Rodriquez rotation formula

end