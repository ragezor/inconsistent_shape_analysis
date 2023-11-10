function nu = heat_flow(nu, F2V,V2F,Lap,iter, beta, alpha)
        nu = F2V*nu;
        mu = nu;
        nu_last = mu;
        for j = 1:iter
            nu = nu_last + beta*Lap*nu_last - alpha*(nu_last - mu);
            nu_last = nu;
        end
        nu = V2F*nu;
end