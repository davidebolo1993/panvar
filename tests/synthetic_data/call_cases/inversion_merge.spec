# An inversion site with two inversions that start at one breakpoint: imid inverts X2 X3 (1400 bp), imid2 X2 X3 Xs
# (1650 bp). They share most of their bases, so the merge criteria make them one event: one record, with the range
# of their lengths, each haplotype keeping its own stretch in the carrier table.
N F0 200 1601
N R 300 1602
N X1 400 1603
N X2 400 1604
N X3 1000 1605
N Xs 250 1606
N X4 400 1607
N F9 200 1608
P ref#0#chr1:1-3450 F0+,R+,X1+,X2+,X3+,Xs+,X4+,R-,F9+
P r1#1#chr1 F0+,R+,X1+,X2+,X3+,Xs+,X4+,R-,F9+
P imid#1#chr1 F0+,R+,X1+,X3-,X2-,Xs+,X4+,R-,F9+
P imid2#1#chr1 F0+,R+,X1+,Xs-,X3-,X2-,X4+,R-,F9+
