# A block (A, the 60 bp unit R repeated, B) the reference holds twice, with 3 then 2 copies of R. u1 changes the unit
# count in the first block, u2 in the second, u12 in both, l lacks the second block. The repeat is counted over both
# blocks, one record for the site: 5 7 7 4 3 for r1 u1 u2 u12 l.
N F0 200 1301
N A 100 1302
N R 60 1303
N B 100 1304
N M 150 1305
N F9 200 1306
P ref#0#chr1:1-1250 F0+,A+,R+,R+,R+,B+,M+,A+,R+,R+,B+,F9+
P r1#1#chr1 F0+,A+,R+,R+,R+,B+,M+,A+,R+,R+,B+,F9+
P u1#1#chr1 F0+,A+,R+,R+,R+,R+,R+,B+,M+,A+,R+,R+,B+,F9+
P u2#1#chr1 F0+,A+,R+,R+,R+,B+,M+,A+,R+,R+,R+,R+,B+,F9+
P u12#1#chr1 F0+,A+,R+,B+,M+,A+,R+,R+,R+,B+,F9+
P l#1#chr1 F0+,A+,R+,R+,R+,B+,F9+
B A B R
