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



import logging
logger = logging.getLogger(__file__)


import io
import os
import sys
import gzip
import hashlib

from Bio import SeqIO
import Bio

class SeqParser:
    """
    Largely independent module, needs 3 pieces of information passed back in from the outside
    """
    def __init__(self, filepath, num, k):
        if type(filepath) is not str:
            raise TypeError("kdb.seqrecord.SeqParser expects a str as its first positional argument")
        elif type(num) is not int:
            raise TypeError("kmerdb.seqparser.SeqParser expects an int as its second positional argument")
        elif type(k) is not int:
            raise TypeError("kmerdb.seqparser.SeqParser expects an int as its third positional argument")
        
        self.k = k
        self.num = num
        self.reads = []
        # Header items
        self.filepath = filepath
        self.md5 = None
        self.sha256 = None
        self.total_reads = 0
        self.total_kmers = 0
        self.unique_kmers = 0
        self.nullomers = 0
        self.mononucleotides = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        # 
        self.compressed = False
        self.fastq = False
        exts = os.path.splitext(filepath)

        if exts[-1] == ".gz":
            self.compressed = True
            nogzexts = os.path.splitext(exts[0])
            if nogzexts[-1] == ".fq" or nogzexts[-1] == ".fastq":
                self.fastq = True
            elif nogzexts[-1] == ".fna" or nogzexts[-1] == ".fasta" or exts[-1] == ".fa":
                self.fastq = False
            else:
                raise ValueError("Cannot parse files of extension '{0}'.\n\nRequires fasta (.fna, .fasta, .fa), fastq (.fq, .fastq), or their gzipped equivalents")
        else: # Must be fasta or fastq uncompressed
            if exts[-1] == ".fq" or exts[-1] == ".fastq":
                self.fastq = True
            elif exts[-1] == ".fna" or exts[-1] == ".fasta" or exts[-1] == ".fa":
                self.fastq = False
            else:
                raise ValueError("Cannot parse files of extension '{0}'.\n\nRequires fasta (.fna, .fasta, .fa), fastq (.fq, .fastq), or their gzipped equivalents")
        if self.fastq:
            self.__class__.__iter__ = self._iter_fastq
            self.__class__.__next__ = self._next_fastq
        else:
            self.__class__.__iter__ = self._iter_fasta
            self.__class__.__next__ = self._next_fasta

                

                
        if self.compressed:
            if self.fastq:
                logger.info("Opening gzipped fastq file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fastq")
            else:
                logger.info("Opening gzipped fasta file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(gzip.open(self.filepath, 'rt'), "fasta")
        else:
            if self.fastq:
                logger.info("Opening uncompressed fastq file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(open(self.filepath, 'r'), "fastq")
            else:
                logger.info("Opening uncompressed fasta file '{0}'...".format(filepath))
                self._handle = SeqIO.parse(open(self.filepath, 'r'), "fasta")

        # Get checksums
        self.md5, self.sha256 = self.__checksum()


    def __checksum(self):
        """Generates md5 and sha256 checksums of a file
        :returns: (md5, sha256)
        :rtype: tuple
        """
        if not os.path.exists(self.filepath):
            raise IOError("kmerdb.seqparser.FastqParser.__checksum could not find '{}' on the filesystem".format(filepath))
        hash_md5 = hashlib.md5()
        hash_sha256 = hashlib.sha256()
        with open(self.filepath, 'rb') as ifile:
            for chunk in iter(lambda: ifile.read(4096), b""):
                hash_md5.update(chunk)
                hash_sha256.update(chunk)
        return (hash_md5.hexdigest(), hash_sha256.hexdigest())

    def header_dict(self):
        """ Create a header dictionary to convert into YAML to go in the header block of the compression header. Has a schema to be validated, defined in config.py

        :returns: dict
        :rtype: 

        """
        return {
            "filename": self.filepath,
            "md5": self.md5,
            "sha256": self.sha256,
            "total_reads": self.total_reads,
            "total_kmers": self.total_kmers,
            "unique_kmers": self.unique_kmers or 4**self.k - self.nullomers,
            "nullomers": self.nullomers,
            "mononucleotides": self.mononucleotides
        }
            
    def __exit__(self, exc_type, exc_value, traceback):
        self._handle.close()

    def __enter__(self):
        return self
        
        
    def _iter_fastq(self):
        try:
            for i in range(self.num):
                self.reads.append(self._handle.__next__())
        except StopIteration as e:
            # with db._engine.connect() as conn:
            #     self.nullomers = db.get_nullomers(conn).pop()
            #     self.total_kmers = db.get_sum_counts(conn).pop()
            #     self.unique_kmers = 4**k - self.nullomers
            pass
        except ValueError as e:
            logger.error(e)
            logger.error("\n\nFastq format error: '{0}' seems to not be fastq format\n\n")
            raise e
        return self

    def _next_fastq(self):
        if not len(self.reads):
            raise StopIteration
        else:
            self.total_reads += 1
            #sys.stderr.write("\r")
            #sys.stderr.write("Read {0} reads from '{1}'...".format(self.total_reads, self.filepath))
            read = self.reads.pop()
            self.mononucleotides["A"] += read.seq.count('A')
            self.mononucleotides["T"] += read.seq.count('T')
            self.mononucleotides["C"] += read.seq.count('C')
            self.mononucleotides["G"] += read.seq.count('G')
            # If you needed to see the fastq mononucleotides go up during parsing
            #sys.stderr.write("\r")
            #sys.stderr.write("Mononucleotides | A: {0}, C: {1}, G: {2}, T{3}".format(self.mononucleotides["A"], self.mononucleotides["C"], self.mononucleotides["G"], self.mononucleotides["T"]))
            return read

    def _iter_fasta(self):
        return self

    def _next_fasta(self):

        seq = next(self._handle)
        self.mononucleotides["A"] += seq.seq.count('A')
        self.mononucleotides["T"] += seq.seq.count('T')
        self.mononucleotides["C"] += seq.seq.count('C')
        self.mononucleotides["G"] += seq.seq.count('G')

        logger.debug("Mononucleotides | A: {0}, C: {1}, G: {2}, T: {3}".format(self.mononucleotides["A"], self.mononucleotides["C"], self.mononucleotides["G"], self.mononucleotides["T"]))
        self.total_reads += 1
        sys.stderr.write("\r")
        sys.stderr.write("Read {0} sequences from '{1}'...\n".format(self.total_reads, self.filepath))
        
        return seq
        
        

