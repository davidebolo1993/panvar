# An inversion of M1 to M5 (Q among them). inv carries it cleanly; invq walks Q, a short node of the inverted
# stretch, forwards in its left flank instead, so that flank reaches into the stretch. Only Q is left out when the
# runs are tested to lie apart: both haplotypes carry one inversion, at the same breakpoints, although short nodes
# lie at every end of the three runs.
N F0 2000 1701
N A1 1000 1702
N A2 40 1703
N M1 40 1704
N M2 600 1705
N Q 60 1706
N M3 600 1707
N M4 600 1708
N M5 40 1709
N B1 40 1710
N F9 3000 1711
P ref#0#chr1:1-8020 F0+,A1+,A2+,M1+,M2+,Q+,M3+,M4+,M5+,B1+,F9+
P r1#1#chr1 F0+,A1+,A2+,M1+,M2+,Q+,M3+,M4+,M5+,B1+,F9+
P inv#1#chr1 F0+,A1+,A2+,M5-,M4-,M3-,Q-,M2-,M1-,B1+,F9+
P invq#1#chr1 F0+,A1+,Q+,A2+,M5-,M4-,M3-,M2-,M1-,B1+,F9+
