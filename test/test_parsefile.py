import sys
import os
import unittest

sys.path.append("..")

from kmerdb import parse, database

class Test_parsefile(unittest.TestCase):

    def setUp(self):

        self.k = 3
        self.fasta  = "test/data/Ctetani_E88.fasta.gz"


    def test_str_filepath_only(self):
        """
        kmer_to_id throws a TypeError if it receives anything but a str
        """
        with self.assertRaises(TypeError):
            parse.parsefile(None, self.k)
        with self.assertRaises(TypeError):
            parse.parsefile(1, self.k)
        with self.assertRaises(TypeError):
            parse.parsefile(1.0, self.k)
        with self.assertRaises(TypeError):
            parse.parsefile([1], self.k)
        with self.assertRaises(TypeError):
            parse.parsefile({'hello': 'world'}, self.k)

    def test_filepath_exists(self):
        """
        parsefile throws an OSError if it receives a filepath that doesnt exist
        """
        with self.assertRaises(OSError):
            db, header = parse.parsefile("foobar.txt", self.k)
            db.conn.close()
            os.unlink("foobar.txt")

    def test_k_is_int_only(self):
        """
        parsefile throws a TypeError if it receives anything but an int
        """
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, None)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, True)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, 'hello')
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, 1.0)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, [1])
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, {'hello': 'world'})
        
    def test_p_is_int_only(self):
        """
        parsefile throws a TypeError if it receives anything but an int
        """
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p=None)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p=True)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p='hello')
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p=1.0)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p=[1])
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, p={'hello': 'world'})

    def test_b_is_int_only(self):
        """
        parsefile throws a TypeError if it receives anything but an int
        """
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b=None)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b=True)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b='hello')
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b=1.0)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b=[1])
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, b={'hello': 'world'})

    def test_stranded_is_bool_only(self):
        """
        parsefile throws a TypeError if it receives anything but an int
        """
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded=None)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded='hello')
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded=1)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded=1.0)
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded=[1])
        with self.assertRaises(TypeError):
            parse.parsefile(self.fasta, self.k, stranded={'hello': 'world'})


    def test_returns_database_and_header_dictionary(self):
        """
        parsefile returns a header dictionary object
        """
        testdb, headerDict = parse.parsefile(self.fasta, self.k)
        self.assertDictEqual({'filename': 'test/data/Ctetani_E88.fasta.gz',
                              'md5': '3a8682befcdba8ca6405a4bc238e6f06',
                              'mononucleotides': {'A': 1016563, 'C': 418556, 'G': 404380, 'T': 1033834},
                              'nullomers': 0,
                              'sha256': '438446b1081c6aa65c88cc81fbe5e65c85077f1c46951e52793f8f4bc5fe7b67',
                              'total_kmers': 5746658,
                              'total_reads': 2,
                              'unique_kmers': 64}
                             , headerDict)
        self.assertTrue(isinstance(testdb, database.SqliteKdb))        
            

if __name__ == '__main__':
    unittest.main()
