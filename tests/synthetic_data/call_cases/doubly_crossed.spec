# Sites the reference crosses twice. One copy (A, a site with B or B2, C, a site with E or E2, F) lies twice on the
# reference with a unique stretch (M, N) between, so each site's walk runs from its first crossing to its last and
# takes in the stretch. Haplotypes i* carry an insertion I in the stretch: both sites see it, and it must be
# written once. Haplotypes a* carry the other allele of both sites in the first copy.
N S 200 301
N A 60 302
N B 60 303
N B2 60 304
N C 60 305
N E 60 306
N E2 60 307
N F 60 308
N M 300 309
N I 120 310
N N 300 311
N T 200 312
P ref#0#chr1:1-1600 S+,A+,B+,C+,E+,F+,M+,N+,A+,B+,C+,E+,F+,T+
P i1#1#chr1 S+,A+,B+,C+,E+,F+,M+,I+,N+,A+,B+,C+,E+,F+,T+
P i2#1#chr1 S+,A+,B+,C+,E+,F+,M+,I+,N+,A+,B+,C+,E+,F+,T+
P i3#1#chr1 S+,A+,B+,C+,E+,F+,M+,I+,N+,A+,B+,C+,E+,F+,T+
P a1#1#chr1 S+,A+,B2+,C+,E2+,F+,M+,N+,A+,B+,C+,E+,F+,T+
P a2#1#chr1 S+,A+,B2+,C+,E2+,F+,M+,N+,A+,B+,C+,E+,F+,T+
P r1#1#chr1 S+,A+,B+,C+,E+,F+,M+,N+,A+,B+,C+,E+,F+,T+
B A C B;B2
B C F E;E2
