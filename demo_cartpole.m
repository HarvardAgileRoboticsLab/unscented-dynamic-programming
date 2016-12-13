function demo_cartpole
% A cart pole swing-up demo.
clc;
close all;

dt = .1;% time step

% quadratic costs
Qf = 1000*eye(4);
Q = .1*eye(4);
R = .01;

% optimization problem
T  = 5/dt; % horizon
x0 = [0 0 0 0]';
xg = [0 pi 0 0]';
u0 = zeros(1,T); % initial controls

% run the optimization
Op.parallel = false;
Op.tolFun = 1e-10;
Op.tolGrad = 1e-10;
Op.maxIter = 150;
xhistd = 0;
uhistd = 0;
costd = 0;
traced = struct();

[xhisti,uhisti,~,~,~,costi,tracei] = iLQG(@(x,u,i) cartpole_dyn_cst(x,u), x0, u0, Op);
[xhistd,uhistd,~,~,~,costd,traced] = iLQG(@(x,u,i) cartpole_dyn_cst3(x,u), x0, u0, Op);
[xhistu,uhistu,~,~,~,costu,traceu] = UDP(@(x,u,i) cartpole_dyn_cst2(x,u), x0, u0, Op, @cartpole_dynamics, dt, .01);

% dynamics
function [xn, A, B] = update(x,u)
    %4th order Runge-Kutta Step
    xdot1 = cartpole_dynamics(0,x,u);
    xdot2 = cartpole_dynamics(0,x+.5*dt*xdot1,u);
    xdot3 = cartpole_dynamics(0,x+.5*dt*xdot2,u);
    xdot4 = cartpole_dynamics(0,x+dt*xdot3,u);
    
    xn = x + (dt/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
    
    if nargout > 1
        n = length(x);
        m = length(u);
        
        delta = 1e-7;
        Dx = delta*eye(n);
        Du = delta*eye(m);
        
        A1 = zeros(n,n);
        A2 = A1;
        A3 = A1;
        A4 = A1;
        B1 = zeros(n,m);
        B2 = B1;
        B3 = B1;
        B4 = B1;
        for j = 1:n
            xp1 = cartpole_dynamics(0,x+Dx(:,j),u);
            xm1 = cartpole_dynamics(0,x-Dx(:,j),u);
            A1(:,j) = (xp1-xm1)/(2*delta);
            
            xp2 = cartpole_dynamics(0,x+.5*dt*xdot1+Dx(:,j),u);
            xm2 = cartpole_dynamics(0,x+.5*dt*xdot1-Dx(:,j),u);
            A2(:,j) = (xp2-xm2)/(2*delta);
            
            xp3 = cartpole_dynamics(0,x+.5*dt*xdot2+Dx(:,j),u);
            xm3 = cartpole_dynamics(0,x+.5*dt*xdot2-Dx(:,j),u);
            A3(:,j) = (xp3-xm3)/(2*delta);
            
            xp4 = cartpole_dynamics(0,x+dt*xdot3+Dx(:,j),u);
            xm4 = cartpole_dynamics(0,x+dt*xdot3-Dx(:,j),u);
            A4(:,j) = (xp4-xm4)/(2*delta);
        end
        for j = 1:m
            xp1 = cartpole_dynamics(0,x,u+Du(:,j));
            xm1 = cartpole_dynamics(0,x,u-Du(:,j));
            B1(:,j) = (xp1-xm1)/(2*delta);
            
            xp2 = cartpole_dynamics(0,x+.5*dt*xdot1,u+Du(:,j));
            xm2 = cartpole_dynamics(0,x+.5*dt*xdot1,u-Du(:,j));
            B2(:,j) = (xp2-xm2)/(2*delta);
            
            xp3 = cartpole_dynamics(0,x+.5*dt*xdot2,u+Du(:,j));
            xm3 = cartpole_dynamics(0,x+.5*dt*xdot2,u-Du(:,j));
            B3(:,j) = (xp3-xm3)/(2*delta);
            
            xp4 = cartpole_dynamics(0,x+dt*xdot3,u+Du(:,j));
            xm4 = cartpole_dynamics(0,x+dt*xdot3,u-Du(:,j));
            B4(:,j) = (xp4-xm4)/(2*delta);
        end
        
        A = (eye(n)+(dt/6)*A4)*(eye(n)+(dt/3)*A3)*(eye(n)+(dt/3)*A2)*(eye(n)+(dt/6)*A1);
        B = (dt/6)*B4 + (eye(n)+(dt/6)*A4)*(dt/3)*B3 + (eye(n)+(dt/6)*A4)*(eye(n)+(dt/3)*A3)*(dt/3)*B2 + (eye(n)+(dt/6)*A4)*(eye(n)+(dt/3)*A3)*(eye(n)+(dt/3)*A2)*(dt/6)*B1;
        
    end
end

function [A, B] = grad(x,u)
    n = length(x);
    m = length(u);

    A = zeros(n,n);
    B = zeros(n,m);
    
    delta = 1e-7;
    Dx = delta*eye(n);
    Du = delta*eye(m);
    
    for j = 1:n
        xp = update(x+Dx(:,j), u);
        xm = update(x-Dx(:,j), u);
        A(:,j) = (xp-xm)/(2*delta);
    end
    for j = 1:m
        xp = update(x, u+Du(:,j));
        xm = update(x, u-Du(:,j));
        B(:,j) = (xp-xm)/(2*delta);
    end
end

function [fxx, fxu, fuu] = hessian(x,u)
    n = length(x);
    m = length(u);
    
    delta = 1e-5;
    Dx = delta*eye(n);
    Du = delta*eye(m);
    
    fxx = zeros(n,n,n);
    fxu = zeros(n,n,m);
    fuu = zeros(n,m,m);
    
    for j = 1:n
        Ap = grad(x+Dx(:,j), u);
        Am = grad(x-Dx(:,j), u);
        fxx(:,:,j) = (Ap-Am)/(2*delta);
    end
    
    for j = 1:m
        [Ap, Bp] = grad(x, u+Du(:,j));
        [Am, Bm] = grad(x, u-Du(:,j));
        fxu(:,:,j) = (Ap-Am)/(2*delta);
        fuu(:,:,j) = (Bp-Bm)/(2*delta);
    end
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_cst(x,u)

% for a positive-definite quadratic, no control cost (indicated by the 
% iLQG function using nans), is equivalent to u=0
N = size(x,2);
n = size(x,1);
m = size(u,1);

f = zeros(n,N);
c = 0;

if nargout == 2
    for k = 1:N
        if(isnan(u(:,k)))
            c = c + 0.5*(x(:,k)-xg)'*Qf*(x(:,k)-xg);
        else
            f(:,k) = update(x(:,k),u(:,k));
            c = c + 0.5*(x(:,k)-xg)'*Q*(x(:,k)-xg) + 0.5*u(:,k)'*R*u(:,k);
        end
    end
else
    A = zeros(n,n,N);
    B = zeros(n,m,N);

    for k = 1:N
        if(isnan(u(:,k)))
            c = c + 0.5*(x(:,k)-xg)'*Qf*(x(:,k)-xg);
        else
            [f(:,k), A(:,:,k), B(:,:,k)] = update(x(:,k),u(:,k));
            c = c + 0.5*(x(:,k)-xg)'*Q*(x(:,k)-xg) + 0.5*u(:,k)'*R*u(:,k);
        end
    end

    fx  = A(:,:,1:end-1);
    fu  = B(:,:,1:end-1);
    cx  = Q*[x(1,:)-xg(1); x(2,:)-xg(2); x(3,:)-xg(3); x(4,:)-xg(4)];
    cx(:,end) = Qf*(x(:,end)-xg);
    cu  = R*u;
    cxx = repmat(Q, [1 1 N]);
    cxx(:,:,N) = Qf;
    cxu = repmat(zeros(n,m), [1 1 N]);
    cuu = repmat(R, [1 1 N]);
    [f,c,fxx,fxu,fuu] = deal([]); 
end
end

function [f,c,cx,cu,cxx,cxu,cuu] = cartpole_dyn_cst2(x,u)

% for a positive-definite quadratic, no control cost (indicated by the 
% iLQG function using nans), is equivalent to u=0
N = size(x,2);
n = size(x,1);
m = size(u,1);

f = zeros(n,N);
c = 0;

if nargout == 2
    for k = 1:N
        if(isnan(u(:,k)))
            c = c + 0.5*(x(:,k)-xg)'*Qf*(x(:,k)-xg);
        else
            f(:,k) = update(x(:,k),u(:,k));
            c = c + 0.5*(x(:,k)-xg)'*Q*(x(:,k)-xg) + 0.5*u(:,k)'*R*u(:,k);
        end
    end
else
    cx  = Q*[x(1,:)-xg(1); x(2,:)-xg(2); x(3,:)-xg(3); x(4,:)-xg(4)];
    cx(:,end) = Qf*(x(:,end)-xg);
    cu  = R*u;
    cxx = repmat(Q, [1 1 N]);
    cxx(:,:,N) = Qf;
    cxu = repmat(zeros(n,m), [1 1 N]);
    cuu = repmat(R, [1 1 N]);
    [f,c] = deal([]); 
end
end

function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_cst3(x,u)

% for a positive-definite quadratic, no control cost (indicated by the 
% iLQG function using nans), is equivalent to u=0
N = size(x,2);
n = size(x,1);
m = size(u,1);

f = zeros(n,N);
c = 0;

if nargout == 2
    for k = 1:N
        if(isnan(u(:,k)))
            c = c + 0.5*(x(:,k)-xg)'*Qf*(x(:,k)-xg);
        else
            f(:,k) = update(x(:,k),u(:,k));
            c = c + 0.5*(x(:,k)-xg)'*Q*(x(:,k)-xg) + 0.5*u(:,k)'*R*u(:,k);
        end
    end
else
    A = zeros(n,n,N);
    B = zeros(n,m,N);
    fxx = zeros(n,n,n,N);
    fxu = zeros(n,n,m,N);
    fuu = zeros(n,m,m,N);

    for k = 1:N
        if(isnan(u(:,k)))
            c = c + 0.5*(x(:,k)-xg)'*Qf*(x(:,k)-xg);
        else
            f(:,k) = update(x(:,k),u(:,k));
            [A(:,:,k), B(:,:,k)] = grad(x(:,k),u(:,k));
            [fxx(:,:,:,k), fxu(:,:,:,k), fuu(:,:,:,k)] = hessian(x(:,k),u(:,k));
            c = c + 0.5*(x(:,k)-xg)'*Q*(x(:,k)-xg) + 0.5*u(:,k)'*R*u(:,k);
        end
    end

    fx  = A(:,:,1:end-1);
    fu  = B(:,:,1:end-1);
    cx  = Q*[x(1,:)-xg(1); x(2,:)-xg(2); x(3,:)-xg(3); x(4,:)-xg(4)];
    cx(:,end) = Qf*(x(:,end)-xg);
    cu  = R*u;
    cxx = repmat(Q, [1 1 N]);
    cxx(:,:,N) = Qf;
    cxu = repmat(zeros(n,m), [1 1 N]);
    cuu = repmat(R, [1 1 N]);
end
end

end