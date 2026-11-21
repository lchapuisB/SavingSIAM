function Pi_all = build_Pi_all(P)
    J = P.J; Nz = P.Nz;
    Pi_all = cell(J, 1);
    for m = 1:J
        Pi_all{m} = zeros(Nz, Nz);  
        [~ , PI_temp] = discretize_normal(P.Nz, 0, P.rho, sqrt(P.sig2_eps(m)));
        Pi_all{m} = PI_temp'; 
    end
end