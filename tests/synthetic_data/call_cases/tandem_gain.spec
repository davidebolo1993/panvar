# A block of two adjacent sites (A, X or nothing, B; B, I or nothing, C) the reference holds once. g1 and g2 hold a
# second copy right after the first, with I in it: one insertion of the copy, written once, rather than one per
# site. g2 also lacks X in its first copy, a change of that copy of its own.
N F0 200 501
N A 100 502
N X 80 503
N B 100 504
N I 70 505
N C 100 506
N F9 200 507
P ref#0#chr1:1-780 F0+,A+,X+,B+,C+,F9+
P r1#1#chr1 F0+,A+,X+,B+,C+,F9+
P g1#1#chr1 F0+,A+,X+,B+,C+,A+,X+,B+,I+,C+,F9+
P g2#1#chr1 F0+,A+,B+,C+,A+,X+,B+,I+,C+,F9+
B A B X
B B C I
