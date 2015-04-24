function check_input(p,mu,nu,q,mv,nv,r,mw,nw)
% checks compatibility of input parameters

if (nu+p+1 ~= mu)
  error('U, p and Control points dont match!')
end
if (nv+q+1 ~= mv)
  error('V, q and Control points dont match!')
end
if (nw+r+1 ~= mw)
  error('W, r and Control points dont match!')
end