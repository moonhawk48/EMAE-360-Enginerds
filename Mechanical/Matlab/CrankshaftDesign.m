function [R_min, M_min, K_max] = CrankshaftDesign()
    r = 4.32; % crank radius in cm
    w = 5000/60.0; % crank speed in rps
    
    % Defined Variables
    a = 8.826; % cm axial distance between crankpins midpoint to midpoint
    L = 13.239; % cm axial distance between outer pieces of crankweb
    rr = 2.083; % rod ratio
    r_cw = 4.5; % cm radius of crankshaft web around journal bearings
    t = 1.2; % cm thickness of crankshaft web (which will be the same for counterweights)

    m_cs = 3.332; % kg mass of crankshaft web per piston, no counterweight
    m_cr = 1.247; % kg mass of connecting rod
    m_p = 0.952; % kg mass of piston head (including wrist pin)
    
    r_cs = 2.16; % distance from rotational axis to crankshaft web center of mass
    r_cr = 11.831; % distance from wrist pin to connecting rod center of mass

    rho = 7.85; % g/cm^3 for counterweight (4340 from SW)

    % Derived Variables
    l = rr*r*2; % conrod length = rod ratio * stroke length

    m_cp = m_cs*(r_cs/r) + m_cr*(r_cr/l); % approximate mass acting on crankpin
    m_wp = m_p + m_cr*((l-r_cr)/l); % approximate mass acting on wrist pin

    % Iterated Variables
    Ycw = [1 1 1 1 ; 1 0 0 1]; % counterweighting options, four or two

    R = r_cw:0.1:15; % possible options for radius (in cm)
    h = 0:0.1:10; % possible options for height (in cm)
    del = (pi/2):-(pi/360):0; % reasonable options for angle (in degrees)
    alph = pi:-(pi/180):1; % reasonable options for sweep of R (in degrees)

    % Target Values, initiated to ensure they are changed within loop
    minR = 100;
    minM = 10^5;
    maxK = 0;
    maxKR = 100;

    % Iteration
    for R_i = R
        for a_i = alph
            if (R_i*sin(a_i/2) < r_cw) % avoid inward sloping counterweights
                continue;
            end
            for d_i = del
                for h_i = h
                    for i = 1:2
                        [m_cw,g_cw,R_top,R_side] = Counterweights(R_i,h_i,t,a_i,d_i,r_cw,rho); % get add mass and center of mass for current counterweight option
                        if (m_cw < 0 || R_top > l-r) % if counterweight mass is negative or radius conflicts with piston
                            continue;
                        end

                        K = CrankshaftBalancing(w,r,a,L,m_cp,m_wp,m_cw,g_cw,Ycw(i,:));
                        if ~(K >= 75 && K <= 100) % if not within desired balance range
                            continue;
                        end
    
                        R_cw = max(R_top,R_side); % farthest protrusion of counterweight from rotational axis
                        if(sum(Ycw(i,:)*m_cw) < minM) % lightest total counterweight mass at least 75% balanced
                            M_min = ["Min Mass",R_cw,m_cw,K ; R_i,h_i,a_i,d_i ; Ycw(i,:)];
                            minM = sum(Ycw(i,:)*m_cw);
                        end
                        if(K > maxK || (K == maxK && R_cw < maxKR)) % highest balance of moments, optimizing for size secondarily
                            K_max = ["Max Balance",R_cw,m_cw,K ; R_i,h_i,a_i,d_i ; Ycw(i,:)];
                            maxK = K;
                            maxKR = R_cw;
                        end
                        
                        if(R_cw < minR) % smallest outer radius at least 75% balanced
                            R_min = ["Min Size",R_cw,m_cw,K ; R_i,h_i,a_i,d_i ; Ycw(i,:)];
                            minR = R_cw;
                        end
                    end
                end
            end
        end
    end
end