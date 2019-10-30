import io
import sys
import os
import gzip
import tempfile
import yaml
import hashlib
import array
import random
import functools
import copy
import time

import jsonschema
from Bio import SeqIO, bgzf
import boto3
sys.path.append('..')

import config
from kdb import kmer_converter


import logging
logger = logging.getLogger(__file__)

s3 = boto3.resource('s3')
s3client = boto3.client('s3')
s3prefix = "s3://"


class SeqReader:
    def __init__(self, seqFile: str, k: int, forward_only=False):
        if type(seqFile) is not str:
            raise TypeError("kdb.fileutil.SeqReader expects a sequence file (fasta/fastq) filepath or s3 object reference as its first positional argument")
        elif not os.path.exists(seqFile) and seqFile[0:5] != s3prefix:
            raise TypeError("kdb.fileutil.SeqReader expects a sequence file (fasta/fastq) filepath or s3 object reference as its first positional argument")
        elif type(k) is not int:
            raise TypeError("kdb.fileutil.SeqReader expects an integer as its second positional argument")
        elif type(forward_only) is not bool:
            raise TypeError("kdb.fileutil.SeqReader expects the keyword argument 'forward_only' to be a bool")
        self.forward_only = forward_only
        # Generate null profile

        self.k = k
        self.profile = array.array('H')
        for x in range(4**k):
            self.profile.append(0)
        
        logger.info("{} Mb allocaed for k-mer array...".format(round(self.profile.__sizeof__()/1048576, 2)))
        # Download file if S3 object
        if seqFile[0:5] == s3prefix:
            self.filepath = self.__s3_file_download(seqFile)
            self.filepath_is_temporary = True
        else:
            self.filepath = seqFile
        # Get checksums
        self.md5, self.sha256 = self.__checksum(self.filepath)
        # Initialize other attributes
        self.total_kmers, self.unique_kmers, self.total_reads = (0, 0, 0)
        self.is_fasta, self.is_fastq, self.is_gzipped = (False, False, False)
        
        if (self.filepath[-3:] == ".fa" or self.filepath[-6:] == ".fasta"): # It's a fasta file
            self.is_fasta = True
            self.is_fastq = False
            self.is_gzipped = False
            self.__validate_seqfile(self.filepath, fasta=True)
            self.filehandle = open(self.filepath, 'r')
        elif (self.filepath[-6:] == ".fa.gz" or self.filepath[-9:] == ".fasta.gz"): # It's a gzipped fasta file
            self.is_fasta = True
            self.is_fastq = False
            self.is_gzipped = True
            self.__validate_seqfile(self.filepath, fasta=True)
            self.filehandle = open(self.filepath, 'rb')
        elif (self.filepath[-3:] == ".fq" or self.filepath[-6:] == ".fastq"): # It's a fastq file
            self.is_fasta = False
            self.is_fastq = True
            self.is_gzipped = False
            self.__validate_seqfile(self.filepath, fasta=False)
            self.filehandle = open(self.filepath, 'r')
        elif (self.filepath[-6:] == ".fq.gz" or self.filepath[-9:] == ".fastq.gz"): # It's a gzipped fastq file
            self.is_fasta = False
            self.is_fastq = True
            self.is_gzipped = True
            self.__validate_seqfile(self.filepath, fasta=False)
            self.filehandle = open(self.filepath, 'rb')


        
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()
        if self.filepath_is_temporary:
            for f in glob.glob(self.filepath + "*"):
                logger.debug("    Unlinking sequence file '{}'...".format(f))
                os.unlink(f)
            
    def __s3_file_download(self, seqpath):
        """
        Note: the file will be downloaded into a temporary file that needs to be deleted afterwards
        It will create the temporary file with respect to the TMP bash variable 'export TMP=/some/temporary/location'
        :param seqpath: The s3 identifier of a object. 's3://bucket/example.fasta'
        :type seqpath: str
        :returns: The location of a downloaded gennomic Fasta file
        :rtype: str
        """
        if type(seqpath) is not str:
            raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects a str 'seqpath' as its first positional argument")
        elif seqpath[0:5] != s3prefix:
            raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects a s3 object reference its first positional argument. e.g. 's3://bucket/example.txt'")
        seqpath = seqpath.lstrip(s3prefix)
        pathsegs =seqpath.split('/')
        bucket = pathsegs.pop(0)
        key = pathsegs.join('/')
        if seqpath[-3:] == ".gz":
            suffix = '.' + sepath.split('.')[-2:].join('.')
        else:
            suffix = path.splitext(seqpath)[1]
        filepath = tempfile.NamedTemporaryFile(mode='w+b', suffix=suffix, delete=False)
        logger.info("Downloading '{0}' => '{1}'...".format(seqpath, filepath.name))
        obj = s3.Object(bucket, key)
        obj.download_fileobj(filepath)
        filepath.close()
        return filepath.name

    def __validate_seqfile(self, seqpath, fasta=False):
        """

        
        :param seqpath: 
        :returns: {lengths: {seqid: seqlen, ...}, filepath: 'path/to/example.fasta', temporary: False} Returns a dictionary of sequence objects, the filepath, and whether the 
        :rtype: 
        """
        if type(seqpath) is not str:
            raise TypeError("kdb.fileutil.SeqReader.__validate_seqfile expects a str 'seqpath' as its first positional argument")
        elif type(fasta) is not bool:
            raise TypeError("kdb.fileutil.SeqReader.__validate_seqfile expects the keyword argument 'fasta' to be a bool")
        elif not os.path.exists(seqpath) or not os.access(seqpath, os.R_OK):
            raise os.error("kdb.fileutil.SeqReader.__validate_seqfile expects the file to be accessible on the filesystem")
        result = {"lengths": [], "filepath": seqpath}

        try:
            logger.debug("Validating format of sequence file '{0}'...".format(seqpath))
            if fasta is True and (seqpath[-3:] == ".fa" or seqpath[-6:] == ".fa.gz" or seqpath[-6:] == ".fasta" or seqpath[-9:] == ".fasta.gz"): # It's a fasta file
                logger.info("Parsing '{}' as fasta format...".format(seqpath))
                if seqpath[-3:] == ".gz":
                    logger.debug("Gzip compression detected...")
                    with open(result["filepath"], 'rb') as ifile:
                        with gzip.open(ifile, mode='rt') as data:
                            for record in SeqIO.parse(data, "fasta"):
                               pass
                else: # Uncompressed fasta
                    with open(result["filepath"], 'r') as ifile:
                        for record in SeqIO.parse(ifile, "fasta"):
                            pass
            elif fasta is False and (seqpath[-3:] == ".fq" or seqpath[-6:] == ".fq.gz" or seqpath[-6:] == ".fastq" or seqpath[-9:] == ".fastq.gz"): # It's a fastq file
                logger.info("Parsing '{}' as fastq format...".format(seqpath))
                if seqpath[-3:] == ".gz":
                    logger.debug("Gzip compression detected...")
                    with open(result["filepath"], 'rb') as ifile:
                        with gzip.open(ifile, mode='rt') as data:
                            for record in SeqIO.parse(data, "fastq"):
                                pass
                else: # Uncompressed fasta
                    with open(result["filepath"], 'r') as ifile:
                        for record in SeqIO.parse(ifile, "fastq"):
                            pass
            else:
                raise TypeError("kdb.fileutil.SeqReader.__validate_seqfile could not determine the file format of '{}'...".format(seqpath))
        except Exception as e:
            logger.error(e)
            raise e
            #sys.exit(1)
        return None


    def __checksum(self, fname: str):
        """Generates md5 and sha256 checksums of a file

        :param fname: The filepath to generate checksums from
        :type fname: str
        :returns: (md5, sha256)
        :rtype: tuple
        """
        if type(fname) is not str:
            raise TypeError("kdb.fileutil.SeqReader.__checksum expects a filepath str as its argument")
        elif not os.path.exists(fname):
            raise IOError("kdb.fileutil.SeqReader.__checksum could not find '{}' on the filesystem".format(fname))
        hash_md5 = hashlib.md5()
        hash_sha256 = hashlib.sha256()
        with open(fname, 'rb') as ifile:
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
            "total_kmers": self.total_kmers,
            "total_reads": self.total_reads,
            "unique_kmers": self.unique_kmers
        }
    
    def read_profile(self):
        """ Read the profile into memory from the sequencing file. No return, just mutates the object, stores the profile in the .profile attribute.

        :returns: 
        :rtype: NoneType


        """
        i=0
        if self.is_fastq:
            if self.is_gzipped:
                with gzip.open(self.filehandle, mode='rt') as data:
                    for record in SeqIO.parse(data, "fastq"):
                        for c in range(len(record.seq) - self.k + 1):
                            seq = record.seq[c:(c+self.k)]
                            self.profile[kmer_converter.sequence_to_binary(str(seq))] += 1
                            if self.forward_only is False:
                                self.profile[kmer_converter.sequence_to_binary(str(seq.reverse_complement()))] += 1
                        i+=1
            else:
                for record in SeqIO.parse(self.filehandle, "fastq"):
                    for c in range(len(record.seq) - self.k + 1):
                        seq = record.seq[c:(c+self.k)]
                        self.profile[kmer_converter.sequence_to_binary(str(seq))] += 1
                        if self.forward_only is False:
                            self.profile[kmer_converter.sequence_to_binary(str(seq.reverse_complement()))] += 1

                    i+=1
        elif self.is_fasta:
            if self.is_gzipped:
                with gzip.open(self.filehandle, mode='rt') as data:
                    for record in SeqIO.parse(data, "fasta"):
                        for c in range(len(record.seq) - self.k + 1):
                            seq = record.seq[c:(c+self.k)]
                            self.profile[kmer_converter.sequence_to_binary(str(seq))] += 1
                            if self.forward_only is False:
                                self.profile[kmer_converter.sequence_to_binary(str(seq.reverse_complement()))] += 1
                        i+=1
            else:
                for record in SeqIO.parse(self.filehandle, "fasta"):
                    for c in range(len(record.seq) - self.k + 1):
                        seq = record.seq[c:(c+self.k)]
                        self.profile[kmer_converter.sequence_to_binary(str(seq))] += 1
                        if self.forward_only is False:
                            self.profile[kmer_converter.sequence_to_binary(str(seq.reverse_complement()))] += 1
                    i+=1
        else:
            raise TypeError("Could not determine the type of file to parse:\nis_fasta: {0}\nis_fastq: {1}\nis_gzipped: {2}".format(self.is_fasta, self.is_fastq, self.is_gzipped))
        self.total_kmers = functools.reduce(lambda a,b: a+b, self.profile)
        self.unique_kmers = 4**self.k - self.profile.count(0)
        self.total_reads = i
   
class KDBWriter:
    def __init__(self, fileobj: io.IOBase, header: dict):
        if not isinstance(fileobj, io.IOBase):
            raise TypeError("kdb.fileutil.KDBWriter expects a file object as its first positional argument")
        elif type(header) is not dict:
            raise TypeError("kdb.fileutil.KDBWriter expects a dict as its second positional argument")
        try:
            jsonschema.validate(instance=header, schema=config.header_schema)
            if 'body_start' in header:
                header.pop('body_start')
            self.filehandle = bgzf.BgzfWriter(fileobj=fileobj, compresslevel=9)
            self.filepath = fileobj.name
            # Write header
            self.filehandle._write_block(bgzf._as_bytes(yaml.dump(header)))
            self.filehandle.flush()
        except jsonschema.ValidationError as e:
            logger.debug(e)
            raise TypeError("kdb.fileutil.KDBWriter expects valid YAML header data as its first positional argument")
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()
        
    def write_profile(self, profile: array.array, k):
        """ Write the profile to the file defined during KDBWriter __init__. The profile is not bound to this object, is an array.array, this is just the wrapper for the IO operations.
        :param profile: The k-mer profiling data to be written to the .kdb file.
        :type profile: array.array
        :returns: 
        :rtype: NoneType
        
        """
        if not isinstance(profile, array.array):
            raise TypeError("kdb.fileutil.KDBWriter.write_profile expects an array.array as its first argument")
        elif type(k) is not int:
            raise TypeError("kdb.fileutil.KDBWriter.write_profile expects an int as its second positional argument")
        logger.info("Writing profile to  '{}'...".format(self.filepath))
        # Write body
        i=0

        text=''
        size = 0
        for idx, c in enumerate(profile):
            if c > 0:
                neighbors = filter(lambda x: profile[k] > 0, kmer_converter.get_neighbors(idx, k))
                newline = "{0}\t{1}\t{2}\n".format(str(idx), str(c), ','.join(list(map(str, neighbors)))) # FIXME should include neighbor ids as csv
                size += bgzf._as_bytes(newline).__sizeof__()
                i+=1
                if i%1000 == 0:
                    sys.stderr.write("\r")
                    sys.stderr.write("Wrote {} lines to the file...".format(i))
                #logger.debug("Line:\n{0}Size in bytes: {1}".format(newline, size))
                if size > 65534:
                    self.filehandle._write_block(bgzf._as_bytes(text)) # should this be _write_block?
                    self.filehandle.flush()
                    text = ''
                    size = 0
                text += newline

        if text != '': # Final flush of text buffer
            sys.stderr.write("\r")
            sys.stderr.write("Wrote {} lines to the file...\n".format(i))
            logger.info("Writing final block...")
            self.filehandle._write_block(bgzf._as_bytes(text))
            self.filehandle.flush()


class KDBReader:
    def __init__(self, fileobj: io.IOBase):
        if not isinstance(fileobj, io.IOBase):
            raise TypeError("kdb.fileutil.KDBReader expects a file object as its first positional argument")
        self.filehandle = fileobj
        self.filepath   = self.filehandle.name
        self.profile    = array.array('H')
        basename, ext   = os.path.splitext(self.filepath)
        
        if ext == '':
            raise IOError("kdb.fileutil.KDBReader could not determine the extention of '{0}'".format(self.filepath))
        elif ext != '.bgzf' and ext != ".kdb":
            raise IOError("kdb.file.get_header can only parse .bgzf format '.kdb' files.")
        bgzfrdr = bgzf.BgzfReader(fileobj=self.filehandle)
        self.offsets = []
        # Ensure the pointers are at the start of the file
        bgzfrdr.seek(0) 
        self.filehandle.seek(0)
        for values in bgzf.BgzfBlocks(self.filehandle):
            #logger.debug("Raw start %i, raw length %i, data start %i, data length %i" % values)
            self.offsets.append(values) # raw start, raw length, data start, data length
        if len(self.offsets) == 0:
            raise IOError("kdb.fileutil.KDBReader opened an empty file")
        # Ensure the pointers are at the start of the file
        bgzfrdr.seek(0) 
        self.filehandle.seek(0)
        # Read data-length bytes from the zeroth(header) block
        header=bgzfrdr.read(self.offsets[0][3]) 
        logger.info("Read KDB header:\n{0}".format(header))
        data=yaml.safe_load(header)
        if type(data) is not dict:
            raise TypeError("kdb.fileutil.KDBReader the header of '{}' was not parseable YAML formatted".format(self.filepath))
        else:
            try:
                jsonschema.validate(instance=data, schema=config.header_schema)
                data["body_start"] = self.offsets[0][3]+1
                self.header = data
                self.k = self.header['k']
            except jsonschema.ValidationError as e:
                logger.debug(e)
                logger.error("kdb.fileutil.KDBReader couldn't validate the header YAML")
                raise e
            for x in range(4**self.k):
                self.profile.append(0)

            
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()
        
    def read_profile(self):
        """ Read profile from the .kdb file. Stores attribute in the .profile attribute.
        Does not mutate the KDBReader object, just reads data from the file. If the file header doesn't validate, it shouldn't get to this point.
        
        :returns: 
        :rtype: 

        """
        bgzfrdr = bgzf.BgzfReader(fileobj=self.filehandle)
        logger.info("Reading profile from '{}'...".format(self.filepath))
        offsets = copy.deepcopy(self.offsets)
        # Omit the header
        offsets.pop(0)
        # Ensure the pointers are at the start of the file
        bgzfrdr.seek(0) 
        #self.filehandle.seek(0)

        # Read data-length bytes from the zeroth(header) block
        logger.debug("Reading blocks from '{}'...".format(self.filepath))
        for raw_start, raw_length, data_start, data_length in offsets:
            if data_length > 0:
                #logger.debug("Raw start: {0}, Raw length: {1}, Data start: {2}, Data length: {3}".format(raw_start, raw_length, data_start, data_length))
                offset = bgzf.make_virtual_offset(raw_start, 0)
                bgzfrdr.seek(offset)
                #self.filehandle.seek(data_start)
                data = bgzfrdr.read(data_length)

                for line in data.split("\n"):
                    if line != '':
                        # FIXME
                        # This part should include more data definition.
                        # Type validation, neighbor validation, parameterized on
                        # data types, data mutability constraints, eventually subject to the index.
                        idx, count, neighbors = line.split("\t")
                        self.profile[int(idx)] = int(count)

                    
    def add_profile(self, x: array.array):
        if type(x) is not array.array:
            raise TypeError("kdb.fileutil.KDBReader.add_profile expects an array as its argument")
        elif len(self.profile) != len(x):
            raise TypeError("kdb.fileutil.KDBReader.add_profile expects the arrays to have equal length")
        for i in range(len(x)):
            self.profile[i] += x[i]

