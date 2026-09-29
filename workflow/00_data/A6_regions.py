"""Region list for A6/A7 (hs1): 5-Mb chunks classified as fusion_2a / fusion_2b (<=10 Mb from the fusion site),
end (<=10 Mb from a real chromosome end; acrocentric p-arms excluded) or interior (>=30 Mb from ends and
fusion; every 4th chunk). Output: chrom start length class."""
import sys

F = (113_940_058 + 114_049_496) / 2
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}
C = 5_000_000
k = 0
for line in open(sys.argv[1]):
    chrom, ln = line.split()[:2]
    ln = int(ln)
    if not chrom[3:].isdigit():
        continue
    for s in range(0, ln, C):
        l = min(C, ln - s)
        mid = s + l / 2
        if chrom == "chr2" and abs(mid - F) <= 10e6 + C / 2:
            cls = "fusion_2a" if mid < F else "fusion_2b"
        elif s < 10e6 and chrom not in ACRO:
            cls = "end"
        elif s + l > ln - 10e6:
            cls = "end"
        elif min(mid, ln - mid) >= 30e6 and not (chrom == "chr2" and abs(mid - F) < 30e6):
            k += 1
            if k % 4:
                continue
            cls = "interior"
        else:
            continue
        print(chrom, s, l, cls)
