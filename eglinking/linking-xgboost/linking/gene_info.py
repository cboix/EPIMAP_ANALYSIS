#!/usr/bin/env python3
# Author: Benjamin T. James
import sys
import numpy as np

# To save memory, could convert to SQL that has
# offset->mind_list
# mind_list->in
class gene_info:
    def __init__(self, pos_metadata, neg_metadata, dhs_names, log_file=sys.stderr):
        self.log = log_file
        self.pos_list = pos_metadata.loc[pos_metadata['name'].isin(dhs_names)]
        self.neg_list = neg_metadata.loc[neg_metadata['name'].isin(dhs_names)]

        self.pos_list = self.pos_list.assign(
            offset=self.pos_list['tss'].values - self.pos_list['mid'].values)
        # We want the offset to be the difference
        # between tss and mid, so for positive strand:
        self.pos_list.loc[self.pos_list['strand'] == '+', 'offset'] *= -1
        self.pos_list = self.pos_list.sort_values('offset')

    def lookup(self, offset_left, offset_right):
        """Do want left index to be inclusive, but since
        right offset is also used as the next left,
        cannot also have that as inclusive.
        Note: '&' has high precedence, so use parens"""
        print("Loclist")
        loc_list = self.pos_list.loc[(self.pos_list['offset'] >= offset_left) &
                                     (self.pos_list['offset'] < offset_right),
                                     'mind'].drop_duplicates()
        print("positive set lookup")
        pos = self.pos_list.loc[(self.pos_list['mind'].isin(loc_list)) &
                                (self.pos_list['offset'] >= offset_left) &
                                (self.pos_list['offset'] < offset_right)]
        print("negative set lookup")
        neg = self.neg_list.loc[self.neg_list['mind'].isin(loc_list)]
        return pos, neg
