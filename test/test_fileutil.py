import sys
import os
import unittest
import yaml


from kmerdb import fileutil



sys.path.append("..")

filepath_that_doesnt_exist = os.path.join(os.path.dirname(__file__), "data/foo.bar")
valid_kdb_filepath = os.path.join(os.path.dirname(__file__), "data/test_Cac_ATCC824.8.kdb")
valid_kdb_metadata_str = """
count_dtype: uint64
files:
- avg_read_length: 2066440
  filename: test/data/Cacetobutylicum_ATCC824.fasta.gz
  max_read_length: 0
  md5: 0a0f73e1c8b8285703e29279bafaabef
  min_read_length: 50
  nullomers: 1433
  sha256: f9081291b62ff3387f1ca6ee2484669c849ed1840fdf2dd9dc3a0c93e9e87951
  total_kmers: 4132866
  total_reads: 2
  unique_kmers: 64103
frequencies_dtype: float64
header_offset: 377
k: 8
kmer_coverage: 1826
kmer_coverage_histogram_dtype: uint64
kmer_ids_dtype: uint64
metadata: false
metadata_blocks: 1
profile_dtype: uint64
sorted: false
tags: []
total_kmers: 4132866
unique_kmers: 64103
unique_nullomers: 1433
version: 0.8.15
"""
 


class TestDoesInputFileExist(unittest.TestCase):
    def given_a_filepath_that_does_not_exist(self):
        self.filepath_that_doesnt_exist = "data/foo.bar"
        
        
    def then_fileutil_KDBReader_should_throw_a_ValueError(self):

        self.assertRaises(ValueError, fileutil.open, self.filepath_that_doesnt_exist, mode='r', slurp=False)

    def test_nonexistent_kdb_file_with_fileutil_KDBReader(self):
        """
        Actual UnitTest test
        """
        self.given_a_filepath_that_does_not_exist()
        self.then_fileutil_KDBReader_should_throw_a_ValueError()
        



class TestValidKDBFile(unittest.TestCase):


    def given_a_valid_kdb_file(self):
        """ Given a valid .kdb file """
        self.test_kdb_metadata_str = valid_kdb_metadata_str
        self.test_kdb_metadata = yaml.safe_load(self.test_kdb_metadata_str)
        self.test_kdb_filepath = valid_kdb_filepath
        
    def when_checking_if_kdb_file_is_valid(self):
        self.test_valid_kdb = fileutil.open(os.path.join(os.path.dirname(__file__), self.test_kdb_filepath), mode='r', slurp=True)

    def then_yaml_metadata_header_should_be_valid(self):

        
        self.assertTrue(self.test_kdb_metadata == self.test_valid_kdb.metadata)
        self.assertTrue(self.test_valid_kdb.k == 8)

    def then_kmerid_and_count_array_dimensions_should_be_4_power_8(self):
        self.assertTrue(len(self.test_valid_kdb.kmer_ids) == 65536)
        self.assertTrue(len(self.test_valid_kdb.counts) == 65536)

    def test_valid_kdb_file_with_fileutil_KDBReader(self):
        """
        Actual UnitTest test
        """
        self.given_a_valid_kdb_file()
        self.when_checking_if_kdb_file_is_valid()
        self.then_yaml_metadata_header_should_be_valid()
        self.then_kmerid_and_count_array_dimensions_should_be_4_power_8()

    def given_a_filepath_to_an_invalid_kdb_file(self):
        self.filepath_to_invalid_kdb_file = os.path.join(os.path.dirname(__file__), "data/test_invalid_kdb_file.kdb")


    def then_fileutil_open_should_throw_a_ValueError(self):
        # with fileutil.open(self.invalid_kdb_filepath, 'r', slurp=False) as fh:
        #     self.assertRaises(ValueError)
        self.assertRaises(ValueError, fileutil.open, self.filepath_to_invalid_kdb_file, mode='r', slurp=False)
        

    def test_invalid_kdb_file_with_fileutil_open(self):
        """
        Actual UnitTest test
        """

        self.given_a_filepath_to_an_invalid_kdb_file()
        self.then_fileutil_open_should_throw_a_ValueError()



class TestErrorsRaisedByOpen(unittest.TestCase):


    def given_a_filepath_that_is_not_a_string(self):
        self.nonstring_filepath = False
        
    def then_open_raises_TypeError_from_invalid_filepath(self):
        self.assertRaises(TypeError, fileutil.open, self.nonstring_filepath, mode='r', slurp=False)

    def test_invalid_filepath_type(self):
        self.given_a_filepath_that_is_not_a_string()
        self.then_open_raises_TypeError_from_invalid_filepath()

        
    def given_a_mode_that_is_not_a_string(self):
        self.nonstring_mode = False
        self.valid_kdb_filepath = valid_kdb_filepath
        
    def then_open_raises_TypeError_from_invalid_mode_type(self):
        self.assertRaises(TypeError, fileutil.open, self.valid_kdb_filepath, mode=self.nonstring_mode, slurp=False)

    def test_invalid_mode_string_type(self):
        self.given_a_mode_that_is_not_a_string()
        self.then_open_raises_TypeError_from_invalid_mode_type()
        

    def given_a_mode_that_is_not_recognized(self):
        self.invalid_mode_string = 'z'
        self.valid_kdb_filepath = valid_kdb_filepath

    def then_open_raises_ValueError_from_invalid_mode_string(self):
        self.assertRaises(ValueError, fileutil.open, self.valid_kdb_filepath, mode=self.invalid_mode_string, slurp=False)

    def test_invalid_mode_string(self):

        self.given_a_mode_that_is_not_recognized()
        self.then_open_raises_ValueError_from_invalid_mode_string()

    def given_a_slurp_that_is_not_a_bool(self):
        self.invalid_slurp = None
        self.valid_kdb_filepath = valid_kdb_filepath

    def then_open_raises_TypeError_from_invalid_slurp(self):
        self.assertRaises(TypeError, fileutil.open, self.valid_kdb_filepath, mode='r', slurp=self.invalid_slurp)

    def test_invalid_slurp_type(self):
        self.given_a_slurp_that_is_not_a_bool()

        self.then_open_raises_TypeError_from_invalid_slurp()
        

class TestOpenReturnsReaderOrWriterBasedOnMode(unittest.TestCase):

    def given_a_mode_string_to_read_from_filepath(self):
        self.valid_kdb_filepath = valid_kdb_filepath
        self.read_mode = 'r'

    def then_open_returns_KDBReader(self):
        self.assertTrue(isinstance(fileutil.open(self.valid_kdb_filepath, mode=self.read_mode, slurp=False), fileutil.KDBReader))

    def test_open_returns_KDBReader_on_read_mode(self):
        self.given_a_mode_string_to_read_from_filepath()
        self.then_open_returns_KDBReader()
        
    def given_a_mode_string_to_write_to_filepath(self):
        self.valid_kdb_filepath = valid_kdb_filepath
        self.write_mode = 'w'
        self.filepath_that_doesnt_exist = filepath_that_doesnt_exist
        self.test_kdb_metadata_str = valid_kdb_metadata_str
        self.test_kdb_metadata = yaml.safe_load(self.test_kdb_metadata_str)

        
    def then_open_returns_KDBWriter(self):
        self.assertTrue(isinstance(fileutil.open(filepath_that_doesnt_exist, mode=self.write_mode, metadata=self.test_kdb_metadata), fileutil.KDBWriter))

    def test_open_returns_KDBWriter_on_write_mode(self):
        self.given_a_mode_string_to_write_to_filepath()
        self.then_open_returns_KDBWriter()

    def tearDown(self):
        if os.access(filepath_that_doesnt_exist, os.R_OK):
            os.unlink(filepath_that_doesnt_exist)
        

    
if __name__ == "__main__":
    unittest.main()
    
        
