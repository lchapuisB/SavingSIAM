
function B1hat = compute_B1hat(P)
    g = P.gamma; kappa = P.r / (1 + P.r);
    B1hat = g*kappa*P.y_tilde - 0.5*(kappa*g*P.w_tilde).^2 .* P.sig2_eta;
end
