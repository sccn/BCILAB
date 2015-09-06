function arimacom = ssmhtd(model)

%SSMHTD Hillmer-Tiao decomposition of ARIMA state space models.
%   arimacom = SSMHTD(model)

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

% Verify model type and extract subtype information
arimatype   = 'non-arima';
for i = 1 : length(model.info)
    switch model.info{i}.type
        case 'arima'
            arimatype   = 'arima';
            p           = model.info{i}.p;
            d           = model.info{i}.d;
            q           = model.info{i}.q;
            if p == 0 || d == 0, arimatype = 'invalid-arima'; end
        case 'sarima'
            arimatype   = 'sarima';
            p           = model.info{i}.p;
            d           = model.info{i}.d;
            q           = model.info{i}.q;
            P           = model.info{i}.P;
            D           = model.info{i}.D;
            Q           = model.info{i}.Q;
            s           = model.info{i}.s;
            if D == 0 && (p+P == 0 || d == 0), arimatype = 'invalid-arima'; end
        case 'sumarma'
            arimatype   = 'sumarma';
            p           = model.info{i}.p;
            q           = model.info{i}.q;
            D           = model.info{i}.D;
            s           = model.info{i}.s;
            if p == 0 || D == 0, arimatype = 'invalid-arima'; end
        case 'genair'
            arimatype   = 'genair';
            nparam      = model.info{i}.nparam;
            nfreq       = model.info{i}.nfreq;
            freq        = model.info{i}.freq;
        case 'freqspec'
            arimatype   = 'freqspec';
            p           = model.info{i}.p;
            d           = model.info{i}.d;
            q           = model.info{i}.q;
            P           = model.info{i}.P;
            D           = model.info{i}.D;
            nparam      = model.info{i}.nparam;
            nfreq       = model.info{i}.nfreq;
            freq        = model.info{i}.freq;
    end
end

% Retrieve parameter values
switch arimatype
    case 'arima'
        phi     = repmat(NaN, 1, p);
        theta   = repmat(NaN, 1, q);
        zetavar = NaN;
        for i = 1 : model.w
            switch model.paramname{i}
                case 'phi1'
                    if i+p-1 <= model.w && strcmp(model.paramname{i+p-1}, ['phi' int2str(p)])
                        phi     = model.param(i:i+p-1);
                    end
                case 'theta1'
                    if i+q-1 <= model.w && strcmp(model.paramname{i+q-1}, ['theta' int2str(q)])
                        theta   = model.param(i:i+q-1);
                    end
                case 'zeta var'
                    zetavar     = model.param(i);
            end
        end
        if any(isnan([phi theta zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
        TotalPhi        = -phi;
        TotalTheta      = theta;
        D               = 0;
        s               = 12; % Whatever
    case 'sarima'
        phi     = repmat(NaN, 1, p);
        theta   = repmat(NaN, 1, q);
        Phi     = repmat(NaN, 1, P);
        Theta   = repmat(NaN, 1, Q);
        zetavar = NaN;
        for i = 1 : model.w
            switch model.paramname{i}
                case 'phi1'
                    if i+p-1 <= model.w && strcmp(model.paramname{i+p-1}, ['phi' int2str(p)])
                        phi     = model.param(i:i+p-1);
                    end
                case 'theta1'
                    if i+q-1 <= model.w && strcmp(model.paramname{i+q-1}, ['theta' int2str(q)])
                        theta   = model.param(i:i+q-1);
                    end
                case 'Phi1'
                    if i+P-1 <= model.w && strcmp(model.paramname{i+P-1}, ['Phi' int2str(P)])
                        Phi     = model.param(i:i+P-1);
                    end
                case 'Theta1'
                    if i+Q-1 <= model.w && strcmp(model.paramname{i+Q-1}, ['Theta' int2str(Q)])
                        Theta   = model.param(i:i+Q-1);
                    end
                case 'zeta var'
                    zetavar     = model.param(i);
            end
        end
        if any(isnan([phi theta Phi Theta zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
        TotalPhi        = conv([1 -phi], [1 kron(-Phi, [zeros(1, s-1) 1])]);
        TotalTheta      = conv([1 theta], [1 kron(Theta, [zeros(1, s-1) 1])]);
        TotalPhi        = TotalPhi(2:end);
        TotalTheta      = TotalTheta(2:end);
        d               = d + D;
    case 'sumarma'
        phi     = repmat(NaN, 1, p);
        theta   = repmat(NaN, 1, q);
        zetavar = NaN;
        for i = 1 : model.w
            switch model.paramname{i}
                case 'phi1'
                    if i+p-1 <= model.w && strcmp(model.paramname{i+p-1}, ['phi' int2str(p)])
                        phi     = model.param(i:i+p-1);
                    end
                case 'theta1'
                    if i+q-1 <= model.w && strcmp(model.paramname{i+q-1}, ['theta' int2str(q)])
                        theta   = model.param(i:i+q-1);
                    end
                case 'zeta var'
                    zetavar     = model.param(i);
            end
        end
        if any(isnan([phi theta zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
        TotalPhi        = -phi;
        TotalTheta      = theta;
        d               = 0;
    case 'genair'
        a       = NaN;
        b       = NaN;
        c1      = NaN;
        c2      = NaN;
        zetavar = NaN;
        for i = 1 : model.w
            switch model.paramname{i}
                case 'a',           a       = model.param(i);
                case 'b',           b       = model.param(i);
                case 'c1',          c1      = model.param(i);
                case 'c2',          c2      = model.param(i);
                case 'zeta var',    zetavar = model.param(i);
            end
        end
        if nparam == 4
            if any(isnan([a b c1 c2 zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
        else
            if any(isnan([a c1 c2 zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
            b   = -a*c1;
            a   = a + c1;
        end
        mask                    = true(6, 1);
        mask(freq(1:6-nfreq))   = false;
        B                       = [ones(6, 1) [-2*cos((1:5)*pi/6)'; 1] zeros(6, 1)];
        C                       = zeros(6, 1);
        C(mask)                 = c1;
        C(~mask)                = c2;
        B(:, 2)                 = B(:, 2).*C;
        B(1:5, 3)               = C(1:5).^2;
        A                       = [1 -a -b];
        A                       = conv(A, B(1, :));
        A                       = conv(A, B(2, :));
        A                       = conv(A, B(3, :));
        A                       = conv(A, B(4, :));
        A                       = conv(A, B(5, :));
        A                       = conv(A, B(6, 1:2));
        TotalTheta              = A(2:end);
        d           = 1 + 1;
        D           = 1;
        s           = 12;
        TotalPhi    = {};
    case 'freqspec'
        phi     = repmat(NaN, 1, p);
        Phi     = repmat(NaN, 1, P);
        a       = repmat(NaN, 1, q+nparam-3);
        c1      = NaN;
        c2      = NaN;
        zetavar = NaN;
        for i = 1 : model.w
            switch model.paramname{i}
                case 'phi1'
                    if i+p-1 <= model.w && strcmp(model.paramname{i+p-1}, ['phi' int2str(p)])
                        phi     = model.param(i:i+p-1);
                    end
                case 'Phi1'
                    if i+P-1 <= model.w && strcmp(model.paramname{i+P-1}, ['Phi' int2str(P)])
                        Phi     = model.param(i:i+P-1);
                    end
                case 'a1'
                    if i+q-1 <= model.w && strcmp(model.paramname{i+q-1}, ['a' int2str(q)])
                        a       = model.param(i:i+q-1);
                    end
                case 'c1',          c1      = model.param(i);
                case 'c2',          c2      = model.param(i);
                case 'zeta var',    zetavar = model.param(i);
            end
        end
        if any(isnan([phi Phi a c1 c2 zetavar])), error('ssm:ssmhtd:param', 'Missing parameters.'); end
        TotalPhi        = conv([1 -phi], [1 kron(-Phi, [zeros(1, 11) 1])]);
        TotalPhi        = TotalPhi(2:end);
        if nparam == 4, A = [1 -a];
        else A = conv([1 -a], [1 -c1]); end
        mask                    = true(6, 1);
        mask(freq(1:6-nfreq))   = false;
        C                       = zeros(6, 1);
        C(mask)                 = c1;
        C(~mask)                = c2;
        B                       = [ones(6, 1) [-2*cos((1:5)*pi/6)'; 1] zeros(6, 1)];
        B(:, 2)                 = B(:, 2).*C;
        B(1:5, 3)               = C(1:5).^2;
        A                       = conv(A, B(1, :));
        A                       = conv(A, B(2, :));
        A                       = conv(A, B(3, :));
        A                       = conv(A, B(4, :));
        A                       = conv(A, B(5, :));
        A                       = conv(A, B(6, 1:2));
        TotalTheta              = A(2:end);
        d                       = d + D;
        s                       = 12;
    case 'non-arima'
        error('ssm:ssmhtd:notarima', 'model is not arima type.');
    case 'invalid-arima'
        error('ssm:ssmhtd:invalidarima', 'model is not decomposable.');
end

% Hillmer-Tiao Decomposition
[theta ksivar]  = htd(d, D, s, TotalPhi, TotalTheta, zetavar);
arimacom        = ssm_arimacom(d, D, s, TotalPhi, theta, ksivar);

