function xdot = yak_dynamics(t,x,u)

    %State vector:
    %r = x(1:3); %Lab-frame position vector
    r = x(4:6); %MRP rotation from body to lab frame
    v = x(7:9); %Lab-frame velocity vector
    w = x(10:12); %Body-frame angular velocity

    Q = mrptodcm(r);
    
    %Control input:
    thr = u(1); %Throttle command (0-255 as sent to RC controller)
    ail = u(2); %Aileron command (0-255 as sent to RC controller)
    elev = u(3); %Elevator command (0-255 as sent to RC controller)
    rud = u(4); %Rudder command (0-255 as sent to RC controller)

    %Note that body coordinate frame is:
    % x: points forward out nose
    % y: points out right wing tip
    % z: points down

    % ---------- Input Checks ---------- %
    thr = min(255, max(0, thr));
    ail = min(255, max(0, ail));
    elev = min(255, max(0, elev));
    rud = min(255, max(0, rud));

    % ---------- Model Parameters ---------- %
    p = yak_parameters; %load model parameters

    % ---------- Map Control Inputs to Angles ---------- %
    delta_ail = (ail-p.trim_ail)*p.g_ail;
    delta_elev = (elev-p.trim_elev)*p.g_elev;
    delta_rud = (rud-p.trim_rud)*p.g_rud;

    % ---------- Aerodynamic Forces (body frame) ---------- %
    v_body = Q'*v; %body-frame velocity
    v_rout = v_body + cross(w,[0; p.r_ail; 0]);
    v_lout = v_body + cross(w,[0; -p.r_ail; 0]);
    v_rin = v_body + cross(w,[0; p.l_in; 0]) + propwash(thr);
    v_lin = v_body + cross(w,[0; -p.l_in; 0]) + propwash(thr);
    v_elev = v_body + cross(w,[-p.r_elev; 0; 0]) + propwash(thr);
    v_rud = v_body + cross(w,[-p.r_rud; 0; -p.z_rud]) + propwash(thr);

    % --- Outboard Wing Sections --- %
    a_rout = alpha(v_rout);
    a_lout = alpha(v_lout);
    a_eff_rout = a_rout + p.ep_ail*delta_ail; %effective angle of attack
    a_eff_lout = a_lout - p.ep_ail*delta_ail; %effective angle of attack

    F_rout = -p_dyn(v_rout)*.5*p.S_out*[Cd_wing(a_eff_rout); 0; Cl_wing(a_eff_rout)];
    F_lout = -p_dyn(v_lout)*.5*p.S_out*[Cd_wing(a_eff_lout); 0; Cl_wing(a_eff_lout)];

    F_rout = arotate(a_rout,F_rout); %rotate to body frame
    F_lout = arotate(a_lout,F_lout); %rotate to body frame

    % --- Inboard Wing Sections (Includes Propwash) --- %
    a_rin = alpha(v_rin);
    a_lin = alpha(v_lin);
    a_eff_rin = a_rin + p.ep_ail*delta_ail; %effective angle of attack
    a_eff_lin = a_lin - p.ep_ail*delta_ail; %effective angle of attack

    F_rin = -p_dyn(v_rin)*.5*p.S_in*[Cd_wing(a_eff_rin); 0; Cl_wing(a_eff_rin)];
    F_lin = -p_dyn(v_lin)*.5*p.S_in*[Cd_wing(a_eff_lin); 0; Cl_wing(a_eff_lin)];

    F_rin = arotate(a_rin,F_rin); %rotate to body frame
    F_lin = arotate(a_lin,F_lin); %rotate to body frame

    % --- Elevator --- %
    a_elev = alpha(v_elev);
    a_eff_elev = a_elev + p.ep_elev*delta_elev; %effective angle of attack

    F_elev = -p_dyn(v_elev)*p.S_elev*[Cd_elev(a_eff_elev); 0; Cl_plate(a_eff_elev)];

    F_elev = arotate(a_elev,F_elev); %rotate to body frame

    % --- Rudder --- %
    a_rud = beta(v_rud);
    a_eff_rud = a_rud - p.ep_rud*delta_rud; %effective angle of attack

    F_rud = -p_dyn(v_rud)*p.S_rud*[Cd_rud(a_eff_rud); Cl_plate(a_eff_rud); 0];

    F_rud = brotate(a_rud,F_rud); %rotate to body frame

    % --- Thrust --- %
    if thr > p.trim_thr
        F_thr = [(thr-p.trim_thr)*p.g_thr; 0; 0];
        w_mot = [p.g_mot*thr; 0; 0];
    else %deadband
        F_thr = [0; 0; 0];
        w_mot = [0; 0; 0];
    end

    % ---------- Aerodynamic Torques (body frame) ---------- %

    T_rout = cross([0; p.r_ail; 0],F_rout);
    T_lout = cross([0; -p.r_ail; 0],F_lout);

    T_rin = cross([0; p.l_in; 0],F_rin);
    T_lin = cross([0; -p.l_in; 0],F_lin);

    T_elev = cross([-p.r_elev; 0; 0],F_elev);

    T_rud = cross([-p.r_rud; 0; -p.z_rud],F_rud);

    % ---------- Add Everything Together ---------- %

    F_aero = F_rout + F_lout + F_rin + F_lin + F_elev + F_rud + F_thr;
    F = Q*F_aero - [0; 0; p.m*p.g];

    T = T_rout + T_lout + T_rin + T_lin + T_elev + T_rud + cross((p.J*w + p.Jm*w_mot),w);

    xdot = [v;
            .25*((1-r'*r)*w - 2*cross(w,r) + 2*(w'*r)*r);
            F/p.m;
            p.Jinv*T];
end

function a = alpha(v)
    %Angle of attack
    a = atan2(v(3),v(1));
end

function b = beta(v)
    %Sideslip angle
    b = atan2(v(2),v(1));
end

function R = mrptodcm(mrp)
    %Converts a vector of Modified Rodrigues Parameters to a Rotation Matrix
    mrp2 = mrp'*mrp;
    S = hat(mrp);

    R = eye(3) + (8*S*S + 4*(1-mrp2)*S)/((1+mrp2)^2);
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end

function rrot = arotate(a,r)
    %Rotate by angle of attack
    rrot = [cos(a) 0  -sin(a);
              0    1    0;
            sin(a) 0  cos(a)]*r;
end

function rrot = brotate(b,r)
    %Rotate by sideslip angle
    rrot = [cos(b) -sin(b) 0;
            sin(b)  cos(b) 0;
              0       0    1]*r;
end

function v = propwash(thr)
    %Propwash wind speed (body frame)
    %Fit from anemometer data taken at tail
    %No significant difference between wing/tail measurements

    trim_thr = 24; %control input for zero thrust (deadband)
    
    if thr > trim_thr
        v = [5.568*thr^0.199 - 8.859; 0; 0];
    else %deadband
        v = [0; 0; 0];
    end
end

function pd = p_dyn(v)
    %Dynamic pressure
    
    p = yak_parameters; %load model parameters
    
    pd = .5*p.rho*(v'*v);
end

function cl = Cl_wing(a)
    %Lift coefficient (alpha in radians)
    %3rd order polynomial fit to glide-test data
    %Good to about +/- 20 degrees

    a = min(pi/2, max(-pi/2, a));

    cl = -27.52*a^3 - .6353*a^2 + 6.089*a;
end

function cl = Cl_plate(a)
    %Lift coefficient (alpha in radians)
    %Ideal flat plate model used for wing and rudder

    a = min(pi/2, max(-pi/2, a));

    cl = 2*pi*a;
end

function cd = Cd_wing(a)
    %Drag coefficient (alpha in radians)
    %2nd order polynomial fit to glide-test data
    %Good to about +/- 20 degrees

    a = min(pi/2, max(-pi/2, a));

    cd = 2.08*a^2 + .0612;
end

function cd = Cd_elev(a)
    %Drag coefficient (alpha in radians)
    %Induced drag for a tapered finite wing
    %From Phillips P.55

    a = min(pi/2, max(-pi/2, a));
    
    p = yak_parameters; %load model parameters
    
    cd = (4*pi*a^2)/p.Ra_elev;
end

function cd = Cd_rud(a)
    %Drag coefficient (alpha in radians)
    %Induced drag for a tapered finite wing
    %From Phillips P.55

    a = min(pi/2, max(-pi/2, a));
    
    p = yak_parameters; %load model parameters
    
    cd = (4*pi*a^2)/p.Ra_rud;
end