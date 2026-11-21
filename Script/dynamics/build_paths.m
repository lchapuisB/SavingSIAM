function paths = build_paths(P_2010, P_2100, T, mode)
% mode = 'step'        : choc permanent instantané à t=1
%
% mode = 'linear'      : transition linéaire en T périodes (1/T par période)
%      
% mode = 'exponential' : transition exponentielle e^{-lambda*t}

    keys = {'w','w_tilde','y_tilde','sig2_eps','sig2_eta', 'A', 'q'};

    switch lower(mode)

        case 'step'
            for k = 1:numel(keys)
                x0 = P_2010.(keys{k})(:);
                xT = P_2100.(keys{k})(:);
                mat = repmat(xT, 1, T);
                mat(:, 1) = x0;         
                paths.(keys{k}) = mat;
            end

        case 'linear'
            for k = 1:numel(keys)
                x0 = P_2010.(keys{k})(:);
                xT = P_2100.(keys{k})(:);
                paths.(keys{k}) = x0 * linspace(1,0,T) + xT * linspace(0,1,T);
            end

        case 'exponential'
            lam = -log(0.05) / (2*T/3);
            t_arr = 0:(T-1);
            w_old = exp(-lam * t_arr);
            for k = 1:numel(keys)
                x0 = P_2010.(keys{k})(:);
                xT = P_2100.(keys{k})(:);
                paths.(keys{k}) = x0 * w_old + xT * (1 - w_old);
            end

        otherwise
            error('Mode inconnu : %s', mode);
    end
    
end