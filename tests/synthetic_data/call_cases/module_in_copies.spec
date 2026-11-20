# A block (S, a module P Q once or twice, K) the reference holds twice, a gene on P. h1 holds the module twice in its
# first block, h2 in its second, hl lacks the second block, hg holds a third. Copy numbers come from one walk over
# both blocks: the module's count, and the gene's, are 2 3 3 1 3 for r1 h1 h2 hl hg.
N F0 200 1501
N S 100 1502
N P 300 1503
N Q 200 1504
N K 100 1505
N M 150 1506
N F9 200 1507
P ref#0#chr9:1-1950 F0+,S+,P+,Q+,K+,M+,S+,P+,Q+,K+,F9+
P r1#1#chr9 F0+,S+,P+,Q+,K+,M+,S+,P+,Q+,K+,F9+
P h1#1#chr9 F0+,S+,P+,Q+,P+,Q+,K+,M+,S+,P+,Q+,K+,F9+
P h2#1#chr9 F0+,S+,P+,Q+,K+,M+,S+,P+,Q+,P+,Q+,K+,F9+
P hl#1#chr9 F0+,S+,P+,Q+,K+,F9+
P hg#1#chr9 F0+,S+,P+,Q+,K+,M+,S+,P+,Q+,K+,M+,S+,P+,Q+,K+,F9+
B S K P;Q
