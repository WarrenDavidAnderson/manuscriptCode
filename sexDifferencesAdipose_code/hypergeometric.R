

#############################################
## hypergeometric tests for gene set overlaps
## 162 DEGs
#############################################

# obesity
q = 92
m = 92 + 70
n = 33784 - m
k = 92 + 70
res = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

# visceral
q = 24
m = 24 + 138
n = 33784 - m
k = 24 + 138
res = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

# multi tissue
q = 70
m = 70 + 92
n = 33784 - m
k = 70 + 92
res = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

# orthologs
q = 34
m = 34 + 128
n = 33784 - m
k = 34 + 88
res = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)


#############################################
## hypergeometric tests for gene set overlaps
## FCcor genes, 3-way (12242 genes total)
#############################################

library(SuperExactTest)

# elements of the Venn
match = paste0("match",1:8)
gtex = paste0("gtex",1:2627)
decode = paste0("decode",1:262)
aagmex = paste0("aagmex",1:558)
gtex.decode = paste0("gtexdecode",1:538)
gtex.aagmex = paste0("gtexaagmex",1:428)
decode.aagmex = paste0("decodeaagmex",1:86)

# individual data sets
d1 = c(match,gtex,gtex.decode,gtex.aagmex)
d2 = c(match,decode,gtex.decode,decode.aagmex)
d3 = c(match,aagmex,gtex.aagmex,decode.aagmex)

# statistical test
x = list(d1, d2, d3)
n = 12242
res = MSET(x,n,lower.tail=FALSE,log.p=FALSE)

#############################################
## hypergeometric tests for gene set overlaps
## FCcor genes, pairwise (12242 genes total)
#############################################

# fat%
q = 2556
m = q + 1046
n = 12242 - m
k = q + 2796
phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

# fat growth
q = 500
m = q + 394
n = 12242 - m
k = q + 4852
phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)

# IR
q = 705
m = q + 376
n = 12242 - m
k = q + 4647
phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
