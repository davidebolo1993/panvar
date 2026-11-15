# Repeat sites. Bubble 1: three insertions at one point, no haplotype carrying two (a clean site, tagged;
# sample s1 carries two, one on each haplotype). Bubbles 2-3: an insertion and a deletion 8 bp apart, both
# carried by s5#1 (a mixed site, left untagged). Bubbles 4-5: the same, the insertion on contig s6#1#chrA
# and the deletion on s6#01#chrB, one haplotype written two ways (mixed, as --diploid joins them).
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
N J2 70 211
N K2 8 212
N D2 60 213
N R2 300 214
P ref#0#chr1:0-1436 L+,A+,M+,K+,D+,R+,K2+,D2+,R2+
P s1#1#chr1 L+,A+,I1+,M+,K+,D+,R+,K2+,D2+,R2+
P s1#2#chr1 L+,A+,I2+,M+,K+,D+,R+,K2+,D2+,R2+
P s2#1#chr1 L+,A+,I1+,M+,K+,D+,R+,K2+,D2+,R2+
P s2#2#chr1 L+,A+,I3+,M+,K+,D+,R+,K2+,D2+,R2+
P s3#1#chr1 L+,A+,I2+,M+,K+,D+,R+,K2+,D2+,R2+
P s3#2#chr1 L+,A+,I3+,M+,K+,D+,R+,K2+,D2+,R2+
P s4#1#chr1 L+,A+,M+,J+,K+,D+,R+,K2+,D2+,R2+
P s4#2#chr1 L+,A+,M+,K+,R+,K2+,D2+,R2+
P s5#1#chr1 L+,A+,M+,J+,K+,R+,K2+,D2+,R2+
P s5#2#chr1 L+,A+,M+,K+,D+,R+,K2+,D2+,R2+
P s6#1#chrA L+,A+,M+,K+,D+,R+,J2+,K2+,D2+,R2+
P s6#01#chrB L+,A+,M+,K+,D+,R+,K2+,R2+
B A M I1;I2;I3
B M K J
B K R D
B R K2 J2
B K2 R2 D2
