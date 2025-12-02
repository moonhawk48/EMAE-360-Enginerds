function K = CrankshaftBalancing(w,r,a,L,m_cp,m_wp,m_cw,g_cw,Ycw)
    Tx = (m_cp+m_wp)*r*(a/2)*(w^2) - (m_cp+m_wp)*r*(-a/2)*(w^2); % primary moment without counterweights
    %Ty = m_cp*r*(r/2)*(w^2) - m_cp*r*(-r/2)*(w^2);

    if (isequal(Ycw,[1 1 1 1])) % counterweight mass acts at crankpin
        Tx_cw = ((m_cp+m_wp)*r - 2*m_cw*g_cw)*(a/2)*(w^2) - ((m_cp+m_wp)*r - 2*m_cw*g_cw)*(-a/2)*(w^2); % primary moment with counterweights
        %Ty_cw = (m_cp*r - m_cw*g_cw)*(r/2)*(w^2) - (m_cp*r - m_cw*g_cw)*(-r/2)*(w^2);
    elseif (isequal(Ycw,[1 0 0 1])) % counterweight mass acts outside of crankpins
        Tx_cw = ((m_cp+m_wp)*r*(a/2) - m_cw*g_cw*(L/2))*(w^2) - ((m_cp+m_wp)*r*(-a/2) - m_cw*g_cw*(-L/2))*(w^2); % primary moment with counterweights
    else
        error('Invalid Ycw input');
    end

    K = 100*(Tx-Tx_cw)/Tx;
    %Ky = 100*(Ty-Ty_cw)/Ty;
end