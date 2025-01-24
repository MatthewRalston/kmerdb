'''
   Copyright 2020 Matthew Ralston

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''
import sys
import os

#import yaml
#from collections import OrderedDict

import numpy as np

from kmerdb import kmer, config, util


"""
q * r = X1
for q, i in queries
       r, j in references

q * r = match1

"""


cpdef double smith_waterman_with_minimizers(long[:] query_fasta_ids, long[:] reference_fasta_ids, char[:] query_fasta_seqs, char[:] reference_fasta_seqs, long[:] query_coords, long[:] reference_coords, long[:] query_minimizers, long[:] reference_minimizers, int minmatches):

    # Okay so the product of the ref_len and query_len in bp ... is Z: The effective search space of a single match

    # Okay so the product of len(query_mins) ~ z1 in query bp and len(ref_mins) ~ z2 in reference bp is Y, kmer ACGT, alphabetical k-space searched.

    # Is the upper boundary and rate limiting step of traversing all queried sequence space.

    # The product of len(query_fastas) dimX and len(ref_fastas) dimY is X the sequence match search space

    # 

    matching = []
    threetuple = None #(query_fasta_id, reference_fasta_id, kmer_id, )
    # if query_minimizers == 0 then was not a minimizer sequence 
    # okay so the minimizer arrays have become an index of kmer_ids of 
    for q, i in enumerate(query_minimizers):
        for r, j in enumerate(reference_minimizers):
            #minimizer match loop
            if q == 0 or r == 0:
                pass
            else:
                minimizer_match = (q, r, qc, rc, kmer_id)
            
            if matching[q][r] is True:
                
            if q == 0 and r == 0:
                pass # Q and R uninitialized
            else:
                if q > 0 and r > 0:
                    case 1 inner set (intersection of)
                if q > 0 and r == 0:
                    pass
                if q==0 and r > 0:
                    outer = "match"
            if query_minimizers[i] == 0 and reference_minimizers[i] == 0:
                pass # Does nothing because uninitialized (shouldn't happen)
            elif 
            
                


    return match array





    

