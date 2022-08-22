#!/usr/bin/env python3
# Author: Benjamin T. James
import sys, os
from get_enhancer_regions import get_site

def run(enhancer_list, R):
    for rline in R:
        L = rline.strip().split()
        sample_name = L[0]
        print("sample ", sample_name)
        if not os.path.isdir(sample_name):
            print("dir %s does not exist" % sample_name)
            sys.exit(1)
        for enhancer in enhancer_list:
            enh_bed = os.path.join(sample_name, enhancer + ".bed")
            with open(enh_bed, 'w') as W:
                get_site(sample_name, [enhancer], W)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: %s E1 E2 .." % sys.argv[0])
        sys.exit(1)
    run(sys.argv[1:], sys.stdin)
