function C = get_point_curv(p,i,u,U,CP)

N=basisfunc(i,p,u,U);
SumNw = 0;
for b = 0:p
  SumNw = N(b+1)*CP(i-p+b,4)+SumNw;
end
C(1) = 0;
C(2) = 0;
C(3) = 0;
for b = 0:p 
  C(1) = N(b+1)*CP(i-p+b,4)*CP(i-p+b,1)/SumNw+C(1);
  C(2) = N(b+1)*CP(i-p+b,4)*CP(i-p+b,2)/SumNw+C(2);
  C(3) = N(b+1)*CP(i-p+b,4)*CP(i-p+b,3)/SumNw+C(3);
end

 
