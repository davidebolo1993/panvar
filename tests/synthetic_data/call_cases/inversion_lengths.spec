# An inversion site (F0 to X3, inverted repeat R at both ends of the stretch) with inversions of three lengths
# between the copies of R. ifirst inverts X1 X2; imid inverts X2 X3 and imid2 X2 X3 Xs, alleles of one inversion
# (they share most of their bases, as two events merged into one record do); iall inverts X1 to X4, past the site's last node. Each haplotype is a
# carrier of the record holding its inversion, and no allele (.) at the others it overlaps.
N F0 200 1601
N R 300 1602
N X1 400 1603
N X2 400 1604
N X3 1000 1605
N Xs 60 1606
N X4 400 1607
N F9 200 1608
P ref#0#chr1:1-3260 F0+,R+,X1+,X2+,X3+,Xs+,X4+,R-,F9+
P r1#1#chr1 F0+,R+,X1+,X2+,X3+,Xs+,X4+,R-,F9+
P ifirst#1#chr1 F0+,R+,X2-,X1-,X3+,Xs+,X4+,R-,F9+
P imid#1#chr1 F0+,R+,X1+,X3-,X2-,Xs+,X4+,R-,F9+
P imid2#1#chr1 F0+,R+,X1+,Xs-,X3-,X2-,X4+,R-,F9+
P iall#1#chr1 F0+,R+,X4-,Xs-,X3-,X2-,X1-,R-,F9+
