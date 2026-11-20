# A copy-number module (A3, A2 twice on the reference) between the two crossings of a site the reference crosses
# twice (A, a site with B or B2, C). The first site's walk runs from its first crossing to its last and takes in
# the module, so both sites see what the haplotypes change there. htest deletes the second A2 and X: one copy
# fewer at the anchor, held by a deletion record of its own, so it reads the reference's count on the module
# record. That deletion is written once; the module record must still read REF_CN for htest.
N F0 2000 81
N A 60 302
N B 60 303
N B2 60 304
N C 60 305
N L 500 91
N A3 20 82
N A2 40 83
N S 40 84
N X 60 85
N R 500 86
N T 2000 312
P ref#0#chr1:1-5580 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hfull#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hgain#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P htest#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,R+,A+,B+,C+,T+
P ha#1#chr1 F0+,A+,B2+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr1#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr2#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr3#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr4#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr5#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr6#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr7#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr8#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr9#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
P hr10#1#chr1 F0+,A+,B+,C+,L+,A3+,A2+,S+,A3+,A2+,X+,R+,A+,B+,C+,T+
B A C B;B2
B L R A3;A2;S;X
