# Four copies of a site (A, X or Xp, B) in the reference, with X in the first and third and Xp in the second and
# fourth. p holds two copies with X, so it is matched to the first and third and lacks the second and the fourth:
# two stretches. q holds one copy, and its loss of three spans both, so the three losses are one copy-number
# record, where p has lost two copies.
N F0 200 801
N A 100 802
N X 80 803
N Xp 80 804
N B 100 805
N F9 200 806
P ref#0#chr1:1-1520 F0+,A+,X+,B+,A+,Xp+,B+,A+,X+,B+,A+,Xp+,B+,F9+
P r1#1#chr1 F0+,A+,X+,B+,A+,Xp+,B+,A+,X+,B+,A+,Xp+,B+,F9+
P p#1#chr1 F0+,A+,X+,B+,A+,X+,B+,F9+
P q#1#chr1 F0+,A+,X+,B+,F9+
B A B X;Xp
