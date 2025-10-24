function [R_min, K_max] = Crankshaft()
    r = 5;
    rcw = 3; % can be a range, but then needs to match up with mrot for each rcw option
    mrot = 20; % needs to be evaluated for crankshaft + conrod big end
    rho = 7.85; % g/cm^3 for 4340 from SW
    Ycw = [1 1 1 1]; %indicates there are four counterweights
    a = 4; % can also be adjusted based on thickness really

    R = 0:10; % reasonable options for radius (in cm)
    h = 0:20; % reasonable options for height (in cm)
    t = 0:0.05:1; % reasonable options for thickness (in cm)
    o = 0:80; % reasonable options for angle (in degrees)

    minR = 10^10;
    maxK = 0;
    for R_i = R %look at this attrocious code lol
        for h_i = h
            for t_i = t
                for o_i = o
                    [R_cw,Kf,Km] = Counterweights(r,rcw,mrot,rho,Ycw,a,R_i,h_i,t_i,o_i);
                    if(R_cw < minR && max(Kf,Km)<=100 && min(Kf,Km) >= 90) %smallest outer radius at least 90% balanced
                        R_min = [R_cw,0,Kf,0,Km ; R_i,h_i,t_i,o_i];
                        minR = R_cw;
                    end
                    if(min(Kf,Km) > maxK && max(Kf,Km) <= 100) %highest balance but don't overbalance
                        K_max = [R_cw,0,Kf,0,Km ; R_i,h_i,t_i,o_i];
                        maxK = min(Kf,Km);
                    end
                end
            end
        end
    end
end