# A block of two sites (A, X or nothing, B; B, I or nothing, C) the reference holds twice, the second copy inverted,
# then a unique stretch and a third site (S, Y or nothing, T). Each copy is read on its own: i1's insertion lies in
# the first copy, d2's deletion in the second, which a window read in one direction never holds. l1 holds one copy
# and lacks the stretch from the first copy to F9: one deletion, written once, and no allele at the second copy.
N F0 200 401
N A 100 402
N X 80 403
N B 100 404
N I 70 405
N C 100 406
N M 150 407
N F9 200 408
N S 100 409
N Y 60 410
N T 100 411
P ref#0#chr1:1-1570 F0+,A+,X+,B+,C+,M+,C-,B-,X-,A-,F9+,S+,Y+,T+
P r1#1#chr1 F0+,A+,X+,B+,C+,M+,C-,B-,X-,A-,F9+,S+,Y+,T+
P i1#1#chr1 F0+,A+,X+,B+,I+,C+,M+,C-,B-,X-,A-,F9+,S+,Y+,T+
P d2#1#chr1 F0+,A+,X+,B+,C+,M+,C-,B-,A-,F9+,S+,Y+,T+
P l1#1#chr1 F0+,A+,X+,B+,C+,F9+,S+,Y+,T+
P y1#1#chr1 F0+,A+,X+,B+,C+,M+,C-,B-,X-,A-,F9+,S+,T+
B A B X
B B C I
B S T Y
