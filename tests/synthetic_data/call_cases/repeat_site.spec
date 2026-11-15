# Repeat sites. Bubble 1: three insertions at one point, no haplotype carrying two (a clean site, tagged;
# sample s1 carries two, one on each haplotype). Bubbles 2-3: an insertion and a deletion 8 bp apart, both
# carried by s5#1 (a mixed site, left untagged).
N L 300 201
N A 100 202
N I1 60 203
N I2 75 204
N I3 90 205
N M 300 206
N J 70 207
N K 8 208
N D 60 209
N R 300 210
P ref#0#chr1:0-1068 L+,A+,M+,K+,D+,R+
P s1#1#chr1 L+,A+,I1+,M+,K+,D+,R+
P s1#2#chr1 L+,A+,I2+,M+,K+,D+,R+
P s2#1#chr1 L+,A+,I1+,M+,K+,D+,R+
P s2#2#chr1 L+,A+,I3+,M+,K+,D+,R+
P s3#1#chr1 L+,A+,I2+,M+,K+,D+,R+
P s3#2#chr1 L+,A+,I3+,M+,K+,D+,R+
P s4#1#chr1 L+,A+,M+,J+,K+,D+,R+
P s4#2#chr1 L+,A+,M+,K+,R+
P s5#1#chr1 L+,A+,M+,J+,K+,R+
P s5#2#chr1 L+,A+,M+,K+,D+,R+
B A M I1;I2;I3
B M K J
B K R D
