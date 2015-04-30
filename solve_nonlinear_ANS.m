function d = solve_nonlinear_ANS(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb,toll,nstep)

nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));
ndof = 3*nu*nv*nw;


df = f/nstep;
F = zeros(ndof,1);
R = zeros(ndof,1);
d = zeros(ndof,1);

for n = 1:nstep
    
    F = F+df;
    R = R-df;
    
    iter = 0;
    % 1. Initial stiffness matrix
    K = stiff_mat_ANS(p,U,q,V,r,W,CP,E,nue,ngauss);
    
    while ( norm(R)/norm(F) > toll )
        
        % 2. reduce stiffness matrix and load vector
        Kred = K;
        Rred = R;
        for i = length(rb):-1:1
            Kred(:,rb(i)) = 0;
            Kred(rb(i),:) = 0;
            Kred(rb(i),rb(i)) = 1;
            Rred(rb(i),1) = 0;
        end
        % 3. reduced displacement vector (solve linear system)
        dred = Kred\-Rred;
       
        d = d + dred;        % update solution 
        
        K = stiff_mat_ANS(p,U,q,V,r,W,CP,E,nue,ngauss) + stiff_mat_grd_sig_grd_ANS(p,U,q,V,r,W,CP,d,E,nue,ngauss);
        Rint = rhs_grd_sigma_ANS(p,U,q,V,r,W,CP,d,E,nue,ngauss);
        
        R = Rint - F;
        
        iter = iter +1 ;
        for i = length(rb):-1:1
            F(rb(i),1) = 0;
            R(rb(i),1) = 0;
        end
        
        fprintf('Time step number %d\n', n)
        fprintf('Iteration number %d\n', iter)
        fprintf('Normalized rerisual %.5e\n', norm(R)/norm(F))
            
    end
    % New control point
    CPd = CPresult(CP,d);
    CP = CPd; % Aggiorno configurazione
    
end

end