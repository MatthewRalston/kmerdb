import sys
import unittest

from Bio import SeqRecord, Seq

sys.path.append("..")

from kmerdb import kmer

class Test_kmer_to_id(unittest.TestCase):

    def test_binary_encoding(self):
        """
        Test that it can convert a k-mer into a 2-bit encoded binary index.
        """
        kmers = ["AA",
                 "AC",
                 "AG",
                 "AT",
                 "CA",
                 "CC",
                 "CG",
                 "CT",
                 "GA",
                 "GC",
                 "GG",
                 "GT",
                 "TA",
                 "TC",
                 "TG",
                 "TT"
        ]
        k_ids = list(range(16))
        
        for x,y in zip(map(kmer.kmer_to_id, kmers), k_ids):
            self.assertEqual(x,y)

    def test_none_ambiguous(self):
        """
        kmer_to_id returns None if an ambiguous base is detected in the fastq/fasta
        """
        self.assertEqual(None, kmer.kmer_to_id("ATCNATC"))

    def test_str_only(self):
        """
        kmer_to_id throws a TypeError if it receives anything but a str
        """
        with self.assertRaises(TypeError):
            kmer.kmer_to_id(None)
        with self.assertRaises(TypeError):
            kmer.kmer_to_id(1)
        with self.assertRaises(TypeError):
            kmer.kmer_to_id(1.0)
        with self.assertRaises(TypeError):
            kmer.kmer_to_id([1])
        with self.assertRaises(TypeError):
            kmer.kmer_to_id({'hello': 'world'})

class Test_Kmers(unittest.TestCase):
    def setUp(self):

        self.kmer = kmer.Kmers(3)
        self.seq  = "AAA"

    def test_init_raises_typeerror(self):
        """
        kmer.Kmers.__init__() raises a TypeError on non-integers
        """
        with self.assertRaises(TypeError):
            kmer.Kmers(None)
        with self.assertRaises(TypeError):
            kmer.Kmers('hello world')
        with self.assertRaises(TypeError):
            kmer.Kmers(1.0)
        with self.assertRaises(TypeError):
            kmer.Kmers([2])
        with self.assertRaises(TypeError):
            kmer.Kmers({'hello': 'world'})

    def test_init_sets_k(self):
        testKmers = kmer.Kmers(3)
        self.assertIs(type(testKmers.k), int)
        self.assertEqual(testKmers.k, 3)
            
    def test_shred(self):
        seqrec = SeqRecord.SeqRecord(Seq.Seq(self.seq, Seq.IUPAC.unambiguous_dna))
        shredded = self.kmer.shred(seqrec)
        self.assertIs(type(shredded), dict)
        

if __name__ == '__main__':
    unittest.main()
