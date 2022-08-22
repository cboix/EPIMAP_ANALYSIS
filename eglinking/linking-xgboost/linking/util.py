#!/usr/bin/env python3
# Author: Benjamin T. James
from collections import deque
import argparse as ap
import h5py
import os
import pandas as pd
import gzip
import pickle

class fixed_queue:
	def __init__(self, capacity):
		self.data = deque([None for x in range(capacity)])

	def add(self, item):
		self.data.pop()
		self.data.appendleft(item)

	def get(self):
		return list(self.data)

def file_type(value):
	if not os.path.isfile(value):
		raise ap.ArgumentTypeError("%s is not a valid file" % value)
	return value

def h5_type(value):
	value = file_type(value)
	try:
		handle = h5py.File(value, 'r')
		handle.close()
	except:
		raise ap.ArgumentTypeError("%s is not a valid HDF5 file" % value)
	return value

def penalty(value):
	if value.lower() in ["l1", "l2"]:
		return value.lower()
	elif value.lower() == "none":
		return None
	else:
		raise ap.ArgumentTypeError("%s is not a valid regression penalty" % value)


def posint(value):
	ival = int(value)
	if ival <= 0:
		raise ap.ArgumentTypeError("%s is not a valid positive integer" % value)
	return ival


def df(value):
	return pd.read_csv(value, sep="\t", header=None).values


def nonnegint(value):
	ival = int(value)
	if ival < 0:
		raise ap.ArgumentTypeError("%s is not a valid non-negative integer" % value)
	return ival

def pickle_gz(value):
	with gzip.open(value, 'rb') as R:
		data = pickle.load(R)
	# except:
 	# 	raise e
#		raise ap.ArgumentTypeError("%s is not a valid pickle.gz file" % value)
	return data

def list_gz(value):
	try:
		data = pd.read_csv(value, sep="\t", header=None).values[:,0]
	except:
		raise ap.ArgumentTypeError("%s is not a valid 1-column tsv file" % value)
	return data
