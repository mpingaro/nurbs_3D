clear;

ts = [0.01 0.05 0.1 0.2];
ps = [5];
rs = [1 2 3];
hs = [32];
lasds = [1 2];
ltype = ['al';'bl'];

for ri = 1:length(rs)
  r = rs(ri);
  for li = 1:length(lasds)
    lasd = lasds(li);
    for ti = 1:length(ts)
      t = ts(ti);
      for pi = 1:length(ps)
        p = ps(pi);
        for hi = 1:length(hs)
          h = hs(hi);
          resS = Example_study1_func(p,r,h,t,lasd);
          name1 = strcat('res/study1Sr',num2str(r),ltype(li,:),'sl',num2str(1/t),'p',num2str(p),'r',num2str(h));
          save (strcat(name1,'.mat'),'resS');
        end
      end
    end
  end
end
