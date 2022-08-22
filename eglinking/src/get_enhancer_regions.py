#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import requests
import gzip

def get_site(sample_name, state_list, out):
    site = "https://personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg19/CALLS/%s_18_CALLS_segments.bed.gz" % (sample_name.upper())
    r = requests.get(site)
    rdata = gzip.decompress(r.content).decode('utf-8')
    for line in rdata.split('\n'):
        L = line.split()
        if len(L) > 3 and L[3] in state_list:
            print("\t".join(L[:3]), file=out)
    return 0

if __name__ == "__main__":
    if len(sys.argv[1:]) <= 1:
        print("Usage: %s enhancer_list enhancer_identifier_1 enhancer_identifier_2 ..." % sys.argv[0])
        sys.exit(1)
    get_site(sys.argv[1], sys.argv[2:], sys.stdout)
