function E = M2E(M, e)

% Convert the mean anomaly to eccentric anomaly. e is
% eccentricity. 

E_new = M + e*sin(M) + ((e^2)/2)*sin(2*M);
err = E_new - e*sin(E_new) - M;

while(abs(err) > 1e-3)
    err = E_new - e*sin(E_new) - M;
    E_last = E_new;
    E_deriv = 1 - e*cos(E_new);
    E_new = E_last - (err/E_deriv);
end

E = E_new;