#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
	"""Check if path is an existing file.
	  :Parameters:
		  path: Path to the file
	"""
	if not os.path.isfile(path):
		if os.path.isdir(path):
			msg = "{0} is a directory".format(path)
		else:
			msg = "{0} does not exist.".format(path)
		raise argparse.ArgumentTypeError(msg)
	return path


def get_arguments():
	"""Retrieves the arguments of the program.
	  Returns: An object that contains the arguments
	"""
	# Parsing arguments
	parser = argparse.ArgumentParser(description=__doc__, usage=
									 "{0} -h"
									 .format(sys.argv[0]))
	parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
						help="Amplicon is a compressed fasta file (.fasta.gz)")
	parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
						help="Minimum sequence length for dereplication")
	parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
						help="Minimum count for dereplication")
	parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
						help="Chunk size for dereplication")
	parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
						help="kmer size for dereplication")
	parser.add_argument('-o', '-output_file', dest='output_file', type=str,
						default="OTU.fasta", help="Output file")
	return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
	"""Lit le fichier <fasta>.fasta passé en paramètre et renvoie toutes les séquence du fichier de taille minimal
	Parametres
	---------
	amplicon_file : fichier fasta
	minseqlen : taille minimal des séquences
	------
	"""
	with open(amplicon_file, "r") as filin:
		seq = ""
		for line in filin:
			if line.startswith('>'):
				if len(seq) >= minseqlen:
					yield seq
				seq = ""
			else:
				seq = seq + line[:-1]
		if len(seq) >= minseqlen:
			yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
	""" Fonction renvoyant des toutes les séquences présente au moins mincount fois
	et de taille minimal mincount dans le fichier fasta amplicon_file pairées avec
	leur nombre d'occurence.

	Parametres
	---------
	amplicon_file : fichier fasta
	minseqlen : taille minimal des séquences
	mincount : nombre minimal d'occurence des séquences
	------
	"""
	seq_list = []
	for seq in read_fasta(amplicon_file, minseqlen):
		seq_list.append(seq)

	for j in Counter(seq_list).most_common():
		if j[1] > mincount:
			yield j

def get_chunks(sequence, chunk_size):
	""" Fonction séparant une sequence en au moins 4 fragments
	de taille chunk_size

	Parametres
	---------
	sequence : une sequence (list de str)
	chunk_sie : taille des fragments
	---------
	Returns
	---------
	chunk_list : Liste contenant tous les fragments de la sequence
	"""
	if len(sequence) < 4*chunk_size:
		raise ValueError

	chunk_list = []
	i=1
	while(i*chunk_size < len(sequence)):
		chunk_list.append(sequence[i*chunk_size-chunk_size:i*chunk_size])
		i += 1

	return chunk_list



def cut_kmer(sequence, kmer_size):
	""" Fonction générant tout les mots de taille kmer_size dans
	une sequence

	Parametres
	---------
	sequence : une sequence (list de str)
	kmer_size : taille des mots
	---------
	"""
	for i in range(len(sequence)-kmer_size+1):
		yield sequence[i:i+kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
	""" Fonction complétant un dictionnaire kmer_dict pris en entrée
	en ajoutant les kmer de taille kmer_size présent dans la sequence
	au dictionnaire s'il sont absent et ajoutant à la clé la valeur
	id_seq 

	Parametres
	---------
	kmer_dict : dictionnaire de kmer
	sequence : une sequence (list de str)
	id_seq : id de la sequence associée
	kmer_size : taille des mots
	---------
	"""
	for seq in cut_kmer(sequence, kmer_size):
		if seq not in kmer_dict:
			kmer_dict[seq] = [id_seq]
		elif id_seq not in kmer_dict[seq]:
			kmer_dict[seq].append(id_seq)
	return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
	""" Prend un dictionnaire ayant pour clé un index de kmer et
	pour valeur une liste d’identifiant des séquences dont ils proviennent,
	une sequence et longueur de kmer: kmer_size. 
	Renvoie la liste de toutes les séquences possédant des kmer
	en commun

	Parametres
	---------
	kmer_dict : dictionnaire de kmer
	sequence : une sequence (list de str)
	kmer_size : taille des mots
	"""
	homologie = []
	for i in cut_kmer(sequence, kmer_size):
		if i in kmer_dict:
			homologie.append(kmer_dict[i])

	flat_list = [item for sublist in homologie for item in sublist]

	if len(flat_list) == 0:
		raise ValueError

	mates = []
	for i in Counter(flat_list).most_common():
		mates.append(i[0])
	return mates



#==============================================================
# Main program
#==============================================================
def main():
	"""
	Main program function
	"""
	# Get arguments
	args = get_arguments()
	for i in dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount):
		print(i)

if __name__ == '__main__':
	main()