%% THIS SEGMENT TAKEN FROM CODE PROVIDED AT:
% https://www.researchgate.net/publication/359950279_Computer_Code_for_Torsional_Vibration_Analysis

if rpm_max > rpm_input(length(rpm_input))
    P(:,length(P(1,:))+1) = P(:,length(P(1,:)));
    rpm_input(length(rpm_input)+1)=rpm_max;
end

%Finding the angle of PCP and starting all pressure traces with PCP at 0 degrees
for h=1:length(P(1,:))
    Ang=1;
    for k=1:length(P(:,1))
        PCP(h)=max(P(:,h));
        if P(k,h)==PCP(h)
            Ang_PCP(h)=Ang;
        end
        Ang=Ang+inc;
    end
end
for h=1:length(P(1,:))
    shift=cast(Ang_PCP(h)*inc,"int16");
    if shift==length(P(:,1))
        shift=1;
    end
    for k=1:length(P(:,1))
        p_zero(k,h)=P(shift,h);
        shift=shift+1;
        if shift==length(P(:,1))+1
            shift=1;
        end
    end
end

%Interpolation for intermediate engine speeds
%Pressure trace matrix
rpm=(rpm_min:int_rpm:rpm_max);
for h=1:length(rpm_input)
    for k=1:length(rpm)
        if rpm_input(h)==rpm(k)
            p_zero_matrix(:,k)=p_zero(:,h);
        end
    end
end
for h=1:length(P(:,1))
    w=1;
    for k=2:length(rpm)-1
        p_zero_matrix(h,k)=((p_zero(h,w)-p_zero(h,w+1))*(rpm(k)-rpm_input(w+1)))/(rpm_input(w)-rpm_input(w+1))+p_zero(h,w+1);
        if rpm_input(w+1)==rpm(k)
            w=w+1;
        end
    end
end

%Interpolation of the shift angles for intermediate engine speeds
%Angles of PCP for every calculated engine speed
for h=1:length(Ang_PCP)
    if Ang_PCP(h)>700
        Ang_PCP(h)=Ang_PCP(h)-721; %Correction for PCP angles BTDC
    end
    for k=1:length(rpm)
        if rpm_input(h)==rpm(k)
            Ang_PCP_out(k)=Ang_PCP(h);
        end
    end
end
w=1;
for k=2:length(Ang_PCP_out)-1
    Ang_PCP_out(k)=round(((Ang_PCP(w)-Ang_PCP(w+1))*(rpm(k)-rpm_input(w+1)))/(rpm_input(w)-rpm_input(w+1))+Ang_PCP(w+1));
    if rpm_input(w+1)==rpm(k)
        w=w+1;
    end
end

%Shifting back the pressure traces to original position (advancing the angle)
for h=1:length(p_zero_matrix(1,:))
    shift=cast((722-Ang_PCP_out(h))/inc,"int16");
    if Ang_PCP_out(h)==1
        shift=1; %Correction for PCP angles at TDC
    end
    if Ang_PCP_out(h)==0
        Ang_PCP_out(h)=-1; %Correction for PCP angles at TDC
    end
    if Ang_PCP_out(h)<0
        shift=-1*Ang_PCP_out(h); %Correction for PCP angles BTDC
    end
    for k=1:length(P(:,1))
        p_out(k,h)=p_zero_matrix(shift,h);
        shift=shift+1;
        if shift==length(P(:,1))+1
            shift=1;
        end
    end
end

%Final correction for actual engine speeds input (redundancy)
for h=1:length(rpm_input)
    for k=1:length(rpm)
        if rpm_input(h)==rpm(k)
            p_out(:,k)=P(:,h);
        end
    end
end

P = p_out;