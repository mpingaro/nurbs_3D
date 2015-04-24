function [d,fs] = solve_nonlinear(p,U,q,V,r,W,CP,E,nue,ngauss,f,rb,toll,nstep)

nu = length(CP(:,1,1,1));
nv = length(CP(1,:,1,1));
nw = length(CP(1,1,:,1));
ndof = 3*nu*nv*nw;


df = f/nstep;
F = zeros(ndof,1);
R = zeros(ndof,1);
d = zeros(ndof,1);
CP0 = CP;           % Initial configurations Omega 0
CPd = CP0;          % Initial configurations Omega t == Omega t
sig = zeros(6,ngauss(1)*ngauss(2)*ngauss(3),8);

% 1. Initial stiffness matrix
K = stiff_mat(p,U,q,V,r,W,CPd,CP0,E,nue,ngauss);
for n = 1:nstep
    
    F = F+df;
    R = R-df;
    
    iter = 0;
    
    while ( norm(R)/norm(F) > toll )
        

        % 2. reduce stiffness matrix and load vector
        for i = length(rb):-1:1
            K(:,rb(i)) = 0;
            K(rb(i),:) = 0;
            K(rb(i),rb(i)) = 1;
            R(rb(i),1) = 0;
        end
        % 3. reduced displacement vector (solve linear system)
        id = K\-R;
        d = d + id;        % update solution 
        
        % New control point (updated configurations)
        CPd = CPresult(CPd,d);
        
        % Compute residual
        [Rint, sig] = rhs_grd_sigma(p,U,q,V,r,W,CPd,CP0,sig,d,E,nue,ngauss);
        % Compute stiffness matrix
        K = stiff_mat(p,U,q,V,r,W,CPd,CP0,E,nue,ngauss) + stiff_mat_grd_sig_grd(p,U,q,V,r,W,CPd,sig,ngauss);
        
        R = Rint - F;
        
        for i = length(rb):-1:1
            R(rb(i),1) = 0;
        end
        iter = iter +1 ;
        
        fprintf('Time step number %d\n', n)
        fprintf('Iteration number %d\n', iter)
        fprintf('Normalized rerisual %.5e\n', norm(R)/norm(F))
        
        plot2in1(p,q,r,U,V,W,CP0,CPd,rb,f)
        
        keyboard
        
    end
    
end

% 5. calculate complete load vector
fs=K*d;

end