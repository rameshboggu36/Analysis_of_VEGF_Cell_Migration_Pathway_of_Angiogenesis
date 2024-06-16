
function dydt = sen_up(~,y)
dydt = zeros(70,1);
k = ones(22, 1);
k_ = ones(15, 1);
v = ones(23, 1);
km = ones(23, 1);

% Extract variables from the state vector 'y'
V = y(1); R = y(2);
VR = y(3); VR2 = y(4);
NCK = y(5); NCKA = y(6); NCKA_VR2 = y(7);
PAK2 = y(8); PAK2A = y(9);
P38 = y(10); P38A = y(11);
MAPKAPK2 = y(12); MAPKAPK2A = y(13);
HSP27 = y(14); HSP27A = y(15);
SHB = y(16); SHBA = y(17); SHBA_VR2 = y(18);
PRAK = y(19); PRAKA = y(20);
FAK = y(21); FAKA = y(22);
PAXLLIN = y(23); PAXILLINA = y(24);


sh1 = y(25); sh2 = y(26); sh3 = y(27);
sh4 = y(28);sh5 = y(29);sh6 = y(30);
sh7 = y(31);sh8 = y(32);sh9 = y(33);
sh10 = y(34);sh11 = y(35);sh12 = y(36);
sh13 = y(37);sh14 = y(38);sh15 = y(39);
sh16 = y(40);sh17 = y(41);sh18 = y(42);
sh19 = y(43);sh20 = y(44);sh21 = y(45); sh22 = y(46); sh23 = y(47);

sp1 = y(48); sp2 = y(49); sp3 = y(50); 
sp4 = y(51); sp5 = y(52); sp6 = y(53); 
sp7 = y(54); sp8 = y(55); sp9 = y(56); 
sp10 = y(57); sp11 = y(58); sp12 = y(59); 
sp13 = y(60); sp14 = y(61); sp15 = y(62); 
sp16 = y(63); sp17 = y(64); sp18 = y(65); 
sp19 = y(66); sp20 = y(67); sp21 = y(68); sp22 = y(69); 
sp23 = y(70);

r1 = k(1) * y(1) * y(2) - k_(1) * y(3);
r2 = k(2) * (y(3)^2) - k_(2) * y(4);
r3 = k(3) * y(4) * y(5) / (km(3) + y(5));
r4 = k(4) * y(6) * y(4) - k_(4) * y(7);
r5 = k(5) * y(6) / (km(5) + y(6));
r6 = k(6) * y(7) * y(8) / (km(6) + y(8));
r7 = v(7) * y(9) / (km(7) + y(9));
r8 = k(8) * y(9) * y(10) / (km(8) + y(10));
r9 = v(9) * y(11) / (km(9) + y(11));
r10 = k(10) * y(11) * y(12) / (km(10) + y(12));
r11 = v(11) * y(13) / (km(11) + y(13));
r12 = k(12) * y(13) * y(14) / (km(12) + y(14));
r13 = v(13) * y(15) / (km(13) + y(15));
r14 = k(14) * y(4) * y(16) / (km(14) + y(16));
r15 = k(15) * y(17) * y(4) - k_(15) * y(18);
r16 = v(16) * y(17) / (km(16) + y(17));
r17 = k(17) * y(11) * y(19) / (km(17) + y(19));
r18 = v(18) * y(20) / (km(18) + y(20));
r19 = k(19) * y(18) * y(21) / (km(19) + y(21));
r20 = k(20) * y(20) * y(21) / (km(20) + y(21));
r21 = v(21) * y(22) / (km(21) + y(22));
r22 = k(22) * y(22) * y(23) / (km(22) + y(23));
r23 = v(23) * y(24) / (km(23) + y(24));

dydt(1) = -r1;
dydt(2) = -r1;
dydt(3) = r1 - r2;
dydt(4) = r2 - r3 - r4 - r13;
dydt(5) = -r3 + r5;
dydt(6) = r3 - r4 - r5;
dydt(7) = r4 - r6;
dydt(8) = -r6 + r7;
dydt(9) = r6 - r7 - r8;
dydt(10) = -r8 + r9;
dydt(11) = r8 - r9 - r10 - r17;
dydt(12) = -r10 + r11;
dydt(13) = r10 - r11 - r12;
dydt(14) = -r12 + r13;
dydt(15) = r12 - r13;
dydt(16) = -r14 + r16;
dydt(17) = r14 - r15 - r16;
dydt(18) = r15 - r19;
dydt(19) = -r17 + r18;
dydt(20) = r17 - r18 - r20;
dydt(21) = -r19 - r20 + r21;
dydt(22) = r19 + r20 - r21 - r22;
dydt(23) = -r22 + r23;
dydt(24) = r22 - r23;

for i = 25:1:47
    dydt(i) = -v(13)*km(13)*y(i)/((km(13) + y(15))^2);
end 

dydt(36) = y(13)*y(14)/(km(12) + y(14)) - v(13)*km(13)*y(36)/((km(13) + y(15))^2);
dydt(37) = -y(15)/(km(13) + y(15)) - v(13)*km(13)*y(37)/((km(13) + y(15))^2);

for i = 48:1:70
    dydt(i) = -v(23)*km(23)*y(i)/((km(23) + y(24))^2);
end 

dydt(69) = y(22)*y(23)/(km(22) + y(23)) - v(23)*km(23)*y(69)/((km(23) + y(24))^2);
dydt(70) = -y(24)/(km(23) + y(24)) - v(23)*km(23)*y(70)/((km(23) + y(24))^2);
end 
