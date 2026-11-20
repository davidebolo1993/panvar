# A block of two adjacent sites (A, X or nothing, B; B, I or nothing, C) the reference holds twice: its first copy
# has X and I, its second neither. hyb holds one copy, with X and without I: its first site is most like the first
# copy and its second most like the second, and on its walk the two meet at B, where the reference holds the rest
# of both copies. hybr is the same haplotype stored on the other strand. c1 and c2 hold one whole copy each.
N F0 200 1201
N A 100 1202
N X 80 1203
N B 100 1204
N I 70 1205
N C 100 1206
N F9 200 1207
P ref#0#chr1:1-1150 F0+,A+,X+,B+,I+,C+,A+,B+,C+,F9+
P r1#1#chr1 F0+,A+,X+,B+,I+,C+,A+,B+,C+,F9+
P hyb#1#chr1 F0+,A+,X+,B+,C+,F9+
P hybr#1#chr1 F9-,C-,B-,X-,A-,F0-
P c1#1#chr1 F0+,A+,X+,B+,I+,C+,F9+
P c2#1#chr1 F0+,A+,B+,C+,F9+
B A B X
B B C I
