\\ 2, 4, 13, 145, 416, 621, 1856, 4730, 5163, 8113, 18260, 142396, 1650399, 6569927, 19975865, 23773865, 371728346, 517406582, 642281555, 864490611, 1352215149, 1503350552,


s=0;
for(k=1,2*10^9,pr=eulerphi(k);if(s<pr,s+=pr,s-=pr);if(s==0,print1(k,", ")))
