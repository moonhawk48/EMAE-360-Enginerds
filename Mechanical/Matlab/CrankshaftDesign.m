function [R_min, M_min, K_max] = CrankshaftDesign()
    r = 5 + 1/3; % crank radius
    w = 5000/60; % crank speed in rps
    
    % Defined Variables
    a = 4; % cm axial distance between crankpins midpoint to midpoint
    L = 6; % cm axial distance between outer pieces of crankweb
    rr = 1.6; % rod ratio
    r_cw = 3; % cm radius of crankshaft web around journal bearings
    t = 0.5; % cm thickness of crankshaft web (which will be the same for counterweights)

    m_cs = 3.5; % kg mass of crankshaft web per piston, no counterweight
    m_cr = 1.5; % kg mass of connecting rod
    m_p = 1.125; % kg mass of piston head (include wrist pin, should be negligible)
    
    r_cs = 0.95; % distance from rotational axis to crankshaft web center of mass
    r_cr = 7.456; %distance from wrist pin to connecting rod center of mass

    rho = 7.85; % g/cm^3 for 4340 from SW

    % Derived Variables
    l = rr*r*2; % conrod length = rod ratio * 2r

    m_cp = m_cs*(r_cs/r) + m_cr*(r_cr/l); % approximate mass acting on crankpin
    m_wp = m_p + m_cr*((l-r_cr)/l); % approximate mass acting on wrist pin

    % Iterated Variables
    Ycw = [1 1 1 1 ; 1 0 0 1]; % counterweighting options, four or two

    R = 0:50; % reasonable options for radius (in cm)
    h = 0:20; % reasonable options for height (in cm)
    o = 0:80; % reasonable options for angle (in degrees)

    % Target Values
    minR = 10^10;
    minM = 10^10;
    maxK = 0;

    % Iteration
    for R_i = R
        for h_i = h
            for o_i = o
                for i = 1:2
                    [m_cw,g_cw,R_cw] = Counterweights(R_i,h_i,t,a,o_i,r_cw,rho); % get add mass and center of mass for current counterweight option
                    K = CrankshaftBalancing(w,r,a,L,m_cp,m_wp,m_cw,g_cw,Ycw(i,:));
                    
                    if (m_cw < 0 || K > 100)
                        continue;
                    end

                    if(sum(Ycw(i,:)*m_cw) < minM && K >= 90) % lightest total counterweight mass at least 90% balanced
                        M_min = ["Min Mass",R_cw,m_cw,K ; R_i,h_i,t,o_i ; Ycw(i,:)];
                        minM = sum(Ycw(i,:)*m_cw);
                    end
                    if(K > maxK) % highest reciprocal balance but don't overbalance
                        K_max = ["Max Balance",R_cw,m_cw,K ; R_i,h_i,t,o_i ; Ycw(i,:)];
                        maxK = K;
                    end
                    
                    if(R_cw < minR && K >= 90) % smallest outer radius at least 90% balanced
                        R_min = ["Min Size",R_cw,m_cw,K ; R_i,h_i,t,o_i ; Ycw(i,:)];
                        minR = R_cw;
                    end
                end
            end
        end
    end
end