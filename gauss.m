function [GP,GW] = gauss(ng)

if (ng==1)
   GP(1) = 0;
   GW(1) = 2;
elseif (ng==2)
   GP(1) = -.5773502691896259;
   GP(2) = .5773502691896259;
   GW(1) = 1;
   GW(2) = 1;
elseif (ng==3)
   GP(1) = -.7745966692414834;
   GP(2) = 0;
   GP(3) = .7745966692414834;
   GW(1) = .5555555555555556;
   GW(2) = .8888888888888889;
   GW(3) = .5555555555555556;
elseif (ng==4)
   GP(1) = -.86113631159405257524;
   GP(2) = -.33998104358485626481;
   GP(3) = .33998104358485626481;
   GP(4) = .86113631159405257524;
   GW(1) = .34785484513745385736;
   GW(2) = .65214515486254614264;
   GW(3) = .65214515486254614264;
   GW(4) = .34785484513745385736;
elseif (ng==5)
   GP(1) = -.90617984593866399282;
   GP(2) = -.53846931010568309105;
   GP(3) = 0;
   GP(4) = .53846931010568309105;
   GP(5) = .90617984593866399282;
   GW(1) = .23692688505618908749;
   GW(2) = .47862867049936646808;
   GW(3) = .56888888888888888888;
   GW(4) = .47862867049936646808;
   GW(5) = .23692688505618908749;
elseif (ng==6)
   GP(1) = -.9324695142031520;
   GP(2) = -.6612093864662645;
   GP(3) = -.2386191860831969;
   GP(4) = .2386191860831969;
   GP(5) = .6612093864662645;
   GP(6) = .9324695142031520;
   GW(1) = .1713244923791703;
   GW(2) = .3607615730481386;
   GW(3) = .4679139345726911;
   GW(4) = .4679139345726911;
   GW(5) = .3607615730481386;
   GW(6) = .1713244923791703;
elseif (ng==7)
   GP(1) = -.9491079123427585;
   GP(2) = -.7415311855993944;
   GP(3) = -.4058451513773972;
   GP(4) = 0;
   GP(5) = .4058451513773972;
   GP(6) = .7415311855993944;
   GP(7) = .9491079123427585;
   GW(1) = .1294849661688697;
   GW(2) = .2797053914892767;
   GW(3) = .3818300505051189;
   GW(4) = .4179591836734694;
   GW(5) = .3818300505051189;
   GW(6) = .2797053914892767;
   GW(7) = .1294849661688697;
elseif (ng==8)
   GP(1) = -.9602898564975362;
   GP(2) = -.7966664774136267;
   GP(3) = -.5255324099163290;
   GP(4) = -.1834346424956498;
   GP(5) = .1834346424956498;
   GP(6) = .5255324099163290;
   GP(7) = .7966664774136267;
   GP(8) = .9602898564975362;
   GW(1) = .1012285362903763;
   GW(2) = .2223810344533745;
   GW(3) = .3137066458778873;
   GW(4) = .3626837833783620;
   GW(5) = .3626837833783620;
   GW(6) = .3137066458778873;
   GW(7) = .2223810344533745;
   GW(8) = .1012285362903763;
elseif (ng==9)
   GP(1) = -.9681602395076261;
   GP(2) = -.8360311073266358;
   GP(3) = -.6133714327005904;
   GP(4) = -.3242534234038089;
   GP(5) = 0;
   GP(6) = .3242534234038089;
   GP(7) = .6133714327005904;
   GP(8) = .8360311073266358;
   GP(9) = .9681602395076261;
   GW(1) = .0812743883615744;
   GW(2) = .1806481606948574;
   GW(3) = .2606106964029354;
   GW(4) = .3123470770400029;
   GW(5) = .3302393550012598;
   GW(6) = .3123470770400028;
   GW(7) = .2606106964029355;
   GW(8) = .1806481606948574;
   GW(9) = .0812743883615744;
elseif (ng==10)
   GP(1) = -0.97390653;
   GP(2) = -0.86506337;
   GP(3) = -0.67940957;
   GP(4) = -0.43339539;
   GP(5) = -0.14887434;
   GP(6) = 0.14887434;
   GP(7) = 0.43339539;
   GP(8) = 0.67940957;
   GP(9) = 0.86506337;
   GP(10) = 0.97390653;
   GW(1) = 0.06667134;
   GW(2) = 0.14945135;
   GW(3) = 0.21908636;
   GW(4) = 0.26926672;
   GW(5) = 0.29552422;
   GW(6) = 0.29552422;
   GW(7) = 0.26926672;
   GW(8) = 0.21908636;
   GW(9) = 0.14945135;
   GW(10) = 0.06667134;
elseif (ng==15)
   GP(1) = -.9879925180204854;
   GP(2) = -.9372733924007059;
   GP(3) = -.8482065834104272;
   GP(4) = -.7244177313601700;
   GP(5) = -.5709721726085388;
   GP(6) = -.3941513470775634;
   GP(7) = -.2011940939974345;
   GP(8) = 0;
   GP(9) = .2011940939974345;
   GP(10) = .3941513470775634;
   GP(11) = .5709721726085388;
   GP(12) = .7244177313601700;
   GP(13) = .8482065834104272;
   GP(14) = .9372733924007059;
   GP(15) = .9879925180204854;
   GW(1) = .03075324199611807;
   GW(2) = .07036604748811134;
   GW(3) = .1071592204671351;
   GW(4) = .1395706779261761;
   GW(5) = .1662692058169852;
   GW(6) = .1861610000155741;
   GW(7) = .1984314853271374;
   GW(8) = .2025782419255562;
   GW(9) = .1984314853271374;
   GW(10) = .1861610000155741;
   GW(11) = .1662692058169852;
   GW(12) = .1395706779261761;
   GW(13) = .1071592204671351;
   GW(14) = .07036604748811134;
   GW(15) = .03075324199611807;
end

end