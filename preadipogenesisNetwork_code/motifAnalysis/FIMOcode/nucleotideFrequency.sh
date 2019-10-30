
# https://bioinformatics.stackexchange.com/questions/361/what-is-the-fastest-way-to-calculate-the-number-of-unknown-nucleotides-in-fasta

# upper vs lowercase (non-repeat vs repeat regions)
# https://bioinformatics.stackexchange.com/questions/225/uppercase-vs-lowercase-letters-in-reference-genome

awk -F a '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 345823800
awk -F c '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 253024277
awk -F g '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 253137432
awk -F t '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 346530271

awk -F A '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 427456346
awk -F C '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 299623606
awk -F G '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 299552717
awk -F T '!/^>"/ {cnt+=NF-1}END{print cnt}' mm10.fa # 427635170


a = 345823800
c = 253024277
g = 253137432
t = 346530271

A = 427456346
C = 299623606
G = 299552717
T = 427635170

# repeat
a / (a+g+c+t) # 0.2885434
c / (a+g+c+t) # 0.2111147
g / (a+g+c+t) # 0.2112091
t / (a+g+c+t) # 0.2891328

# non
A / (A+C+T+G) # 0.2939323
C / (A+C+T+G) # 0.2060306
G / (A+C+T+G) # 0.2059818
T / (A+C+T+G) # 0.2940553

# both
(A+a) / (A+C+T+G+a+g+c+t) # 0.2914976
(C+c) / (A+C+T+G+a+g+c+t) # 0.2083275
(G+g) / (A+C+T+G+a+g+c+t) # 0.2083435
(T+t) / (A+C+T+G+a+g+c+t) # 0.2918314


0.29 + 0.21 + 0.21 + 0.29



