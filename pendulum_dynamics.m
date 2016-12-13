function xdot = pendulum_dynamics(t,x,u)

m = 1;   % kg
l = .5;  % m
b = 0.1; % kg m^2 /s
lc = .5; % m
I = .25; %m*l^2; % kg*m^2
g = 9.81; % m/s^2

q = x(1);
qd = x(2);

qdd = (u - m*g*lc*sin(q) - b*qd)/I;

xdot = [qd; qdd];

end

