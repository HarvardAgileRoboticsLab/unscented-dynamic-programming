function xdot = cartpole_dynamics(t,x,u)

mc = 10;   % mass of the cart in kg
mp = 1;    % mass of the pole (point mass at the end) in kg
l = 0.5;   % length of the pole in m
g = 9.81;  % gravity m/s^2

q = x(1:2);
qd = x(3:4);

s = sin(q(2)); c = cos(q(2));

H = [mc+mp, mp*l*c; mp*l*c, mp*l^2];
C = [0 -mp*qd(2)*l*s; 0 0];
G = [0; mp*g*l*s];
B = [1; 0];

qdd = -H\(C*qd + G - B*u);

xdot = [qd; qdd];

end

