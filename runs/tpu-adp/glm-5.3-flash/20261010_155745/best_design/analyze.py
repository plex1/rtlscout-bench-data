import re
lines = open('vectors.dat').read().splitlines()
rows = []
for ln in lines:
    ln = ln.strip()
    if not ln or ln.startswith('#'):
        continue
    p = ln.split()
    rows.append((int(p[0]), p[1], int(p[2]), p[3], int(p[4]), p[5]))

# Find first tile: first row with w=1 (zero-tile load at line 5)
i = 0
while rows[i][0] == 0:
    i += 1
print("first load start line (1-based):", i+1, "len:", sum(1 for r in rows if r[0]==1 and rows.index(r)>=0))
# first a_valid after load
j = i
while rows[j][2] == 0:
    j += 1
print("first a_valid (line):", j+1)
# first c_valid
k = j
while rows[k][4] == 0:
    k += 1
print("first c_valid (line):", k+1)
print("en-cycles between first a_valid and first c_valid:", k - j)
# count load cycles before it
print("load cycles in that window:", sum(1 for r in rows[i:j] if r[0]==1))

# Now do a general alignment analysis over the whole file:
# hypothesis H(d0): c_valid appears when (num_en_since_accept == d0 + 14)
# where accept = cycle with a_valid=1, num_en counts en-cycles after it.
# Simulate a systolic model to find what the reference does.
N=8
def sim(d0, tapfix='cur'):
    # arrays: w[i][j], acc[j] (col result at bottom), skew regs etc.
    # Instead of full sim, check the alignment using simple counters:
    out = []
    pend = []  # (emit_k, data)
    return out

# Simpler: check hypothesis on identity tile
# find the identity tile: pattern with 7f0000... loads
for idx in range(len(rows)-8):
    if all(rows[idx+t][0]==1 for t in range(8)):
        wdat = [rows[idx+t][1] for t in range(8)]
        if wdat[0].startswith('7f00000000000000'):
            print("identity-ish load start:", idx+1)
            break
j2 = idx+8
while rows[j2][2]==0:
    j2 += 1
print("identity a start:", j2+1)
k2 = j2
while rows[k2][4]==0:
    k2 += 1
print("identity first c:", k2+1, "en-cycles:", k2-j2)
