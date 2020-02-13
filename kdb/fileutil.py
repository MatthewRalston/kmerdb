import io
import sys
import os
import gzip
import tempfile
import yaml, json
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
from kdb import kmer_utils
from kdb import profile as vector_ops



import logging
logger = logging.getLogger(__file__)

s3 = boto3.resource('s3')
s3client = boto3.client('s3')
s3prefix = "s3://"


class SeqReader:
    def __init__(self, seqFile: str, k: int, forward_only=False, keep_temporary=False):
        """Read k-mers from a fasta or fastq file on local storage or streamed from S3 into a temporary file.
        
        """
        # Type check the arguments
        if type(seqFile) is not str:
            raise TypeError("kdb.fileutil.SeqReader expects a sequence file (fasta/fastq) filepath or s3 object reference as its first positional argument")
        elif not os.path.exists(seqFile) and seqFile[0:5] != s3prefix:
            raise TypeError("kdb.fileutil.SeqReader expects a sequence file (fasta/fastq) filepath or s3 object reference as its first positional argument")
        elif type(k) is not int:
            raise TypeError("kdb.fileutil.SeqReader expects an integer as its second positional argument")
        elif type(forward_only) is not bool:
            raise TypeError("kdb.fileutil.SeqReader expects the keyword argument 'forward_only' to be a bool")
        elif type(keep_temporary) is not bool:
            raise TypeError("kdb.fileutil.SeqReader expects the keyword argument 'keep_temporary' to be a bool")

        # Addditional attribute initialization
        self.forward_only = forward_only
        self.k = k
        self.kmer_utils = kmer_utils.Utils(self.k)
        self.total_kmers, self.unique_kmers, self.total_reads = (0, 0, 0)
        self.is_fasta, self.is_fastq, self.is_gzipped = (False, False, False)        
        # Generate null profile
        self.profile = array.array('L')#array.array('H')
        for x in range(4**k):
            self.profile.append(0)
        
        logger.info("{} Mb allocaed for k-mer array...".format(round(self.profile.__sizeof__()/1048576, 2)))
        # Download file if S3 object
        if seqFile[0:5] == s3prefix:
            self.filepath_is_temporary = (not keep_temporary)
            self.filepath = self.__s3_file_download(seqFile, temporary=self.filepath_is_temporary)
        else:
            self.filepath = seqFile
        # Get checksums
        self.md5, self.sha256 = self.__checksum(self.filepath)
        
        if (self.filepath[-3:] == ".fa" or self.filepath[-6:] == ".fasta"): # It's a fasta file
            self.is_fasta, self.is_fastq, self.is_gzipped = True, False, False
            mode = 'r'
        elif (self.filepath[-6:] == ".fa.gz" or self.filepath[-9:] == ".fasta.gz"): # It's a gzipped fasta file
            self.is_fasta, self.is_fastq, self.is_gzipped = True, False, True
            mode = 'rb'
        elif (self.filepath[-3:] == ".fq" or self.filepath[-6:] == ".fastq"): # It's a fastq file
            self.is_fasta, self.is_fastq, self.is_gzipped = False, True, False
            mode = 'r'
        elif (self.filepath[-6:] == ".fq.gz" or self.filepath[-9:] == ".fastq.gz"): # It's a gzipped fastq file
            self.is_fasta, self.is_fastq, self.is_gzipped = False, True, True
            mode = 'rb'
        else:
            raise TypeError("kdb.fileutil.SeqReader.__validate_seqfile could not determine the file format or compression of '{}'...".format(seqpath))
        self.__validate_seqfile(self.filepath, fasta=self.is_fasta)
        self.filehandle = open(self.filepath, mode)
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()

        for f in glob.glob(self.filepath + "*"):
            if self.filepath_is_temporary is True:
                logger.debug("    Unlinking sequence file '{}'...".format(f))
                os.unlink(f)
            else:
                logger.debug("    Not unlinking filepath '{1}'...".format(f))
                
    def __s3_file_download(self, seqpath, temporary=True):
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
        elif type(temporary) is not bool:
            raise TypeError("kdb.fileutil.SeqReader.__s3_file_download expects the keyword argument temporary to be a bool")
        seqpath = seqpath.lstrip(s3prefix)
        pathsegs =seqpath.split('/')
        bucket = pathsegs.pop(0)
        fname = os.path.basename(seqpath)
        key = '/'.join(pathsegs)
        if seqpath[-3:] == ".gz":
            suffix = '.' + '.'.join(seqpath.split('.')[-2:])
        else:
            suffix = path.splitext(seqpath)[1]
        if temporary is True:
            filepath = tempfile.NamedTemporaryFile(mode='w+b', suffix=suffix, delete=False)
            logger.info("Downloading '{0}' => '{1}'...".format(seqpath, filepath.name))
        else:
            filepath = open(fname, 'w+b')
            logger.info("Downloading '{0}' => '{1}'...".format(seqpath, fname))
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
            if fasta is True: #  F A S T A 
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
            elif fasta is False: # F A S T Q
                logger.info("Parsing '{}' as fastq format...".format(seqpath))
                if seqpath[-3:] == ".gz":
                    logger.debug("Gzip compression detected...")
                    with open(result["filepath"], 'rb') as ifile:
                        with gzip.open(ifile, mode='rt') as data:
                            for record in SeqIO.parse(data, "fastq"):
                                pass
                else: # Uncompressed fastq
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
            "unique_kmers": self.unique_kmers,
            "nullomers": self.nullomers
        }
    
    def _shred_seqio_record(self, record):
        """FIXME! briefly describe function

        :param record: a SeqIO record
        :type record: Bio.SeqIO.SeqRecord
        :param j: A total k-mer count to increment
        :type j: int
        :returns: The index of the k-mer, its reverse complement, the total count of k-mers
        :rtype: (int, int, int)

        """
        indexes = []
        for c in range(len(record.seq) - self.k + 1):
            seq = record.seq[c:(c+self.k)]
            idx1, idx2 = self.kmer_utils.sequence_to_binary(str(seq)), self.kmer_utils.sequence_to_binary(str(seq.reverse_complement()))
            if idx1 is None or idx2 is None: # Can only happen if there's a k-mer with a 'N'
                pass
            else:
                indexes.append(idx1)
                if self.forward_only is False:
                    indexes.append(idx2)
                
        return indexes
    
    def read_profile(self):
        """ Read the profile into memory from the sequencing file. No return, just mutates the object, stores the profile in the .profile attribute.

        :returns: 
        :rtype: NoneType


        """


        
        i=1 # Number of fasta sequences / fastq reads
        try:
            if self.is_fastq:
                if self.is_gzipped:
                    with gzip.open(self.filehandle, mode='rt') as data:
                        for record in SeqIO.parse(data, "fastq"):
                            i+=1
                            indexes = self._shred_seqio_record(record)
                            for idx in indexes:
                                self.profile[idx] += 1
                else:
                    for record in SeqIO.parse(self.filehandle, "fastq"):
                            i+=1
                            indexes = self._shred_seqio_record(record)
                            for idx in indexes:
                                self.profile[idx] += 1
            elif self.is_fasta:
                if self.is_gzipped:
                    with gzip.open(self.filehandle, mode='rt') as data:
                        for record in SeqIO.parse(data, "fasta"):
                            i+=1
                            indexes = self._shred_seqio_record(record)
                            for idx in indexes:
                                self.profile[idx] += 1
                else:
                    for record in SeqIO.parse(self.filehandle, "fasta"):
                        i+=1
                        indexes = self._shred_seqio_record(record)
                        for idx in indexes:
                            self.profile[idx] += 1
            else:
                raise TypeError("Could not determine the type of file to parse:\nis_fasta: {0}\nis_fastq: {1}\nis_gzipped: {2}".format(self.is_fasta, self.is_fastq, self.is_gzipped))
        except OverflowError as e:
            logger.warning("Overflow of the C unsigned short type is common when dealing with RNA-Seq, high-volume datasets, or contamination/overrepresentation issues")
            logger.error(e)

            raise OverflowError("kdb.fileutil.SeqReader.read_profile read more than 65,535 of a single k-mer into the profile")
        self.total_kmers = functools.reduce(lambda a,b: a+b, self.profile)
        self.nullomers = self.profile.count(0)
        self.unique_kmers = 4**self.k - self.nullomers
        self.total_reads = i


class KDB:
    def __init__(self, profile: array.array, header: dict, metadata:list=[], neighbors:list=[], make_neighbors:bool=False):
        if not isinstance(profile, array.array):
            raise TypeError("kdb.fileutil.KDB expects an array.array as its first positional arguments")
        elif type(header) is not dict:
            raise TypeError("kdb.fileutil.KDB expects a header dict as its second positional argument")
        elif type(metadata) is not list:
            raise TypeError("kdb.fileutil.KDB expects the keyword argument metadata to be a list")
        elif type(neighbors) is not list:
            raise TypeError("kdb.fileutil.KDB expects the keyword argument neighbors to be a list")
        elif type(make_neighbors) is not bool:
            raise TypeError("kdb.fileutil.KDB expects the keyword argument make_neighbors to be a bool")
        try:
            # More validations and initialization
            jsonschema.validate(instance=header, schema=config.header_schema)
            # if 'body_start' in header:
            #     header.pop('body_start')
        except jsonschema.ValidationError as e:
            logger.debug(e)
            raise TypeError("kdb.fileutil.KDB expects valid YAML header data as its fourth positional argument")
        self.profile       = profile
        self.num_kmers     = len(self.profile)
        self.has_metadata  = False
        self.has_neighbors = False
        self.header        = header
        self.k             = self.header['k']
        self.metadata      = neighbors
        self.neighbors     = metadata
        num_neighbors      = len(self.neighbors)
        num_metadata       = len(self.metadata)

        if 4**self.header['k'] != self.num_kmers:
            raise TypeError("kdb.fileutil.KDB expects the k-mer profile array to have length equal to 4^{0} = {1} but it was {2}...".format(self.header['k'], 4**self.header['k'], num_kmers))
        if num_neighbors != 0 and num_metadata != 0:
            if num_neighbors != num_metadata or num_neighbors != self.num_kmers:
                raise TypeError("kdb.fileutil.KDB expects a k-mer profile array to have equal length to the metadata and neighbors arrays")
            else:
                logger.debug("KDB instance has metadata, neighbors, and k-mer profile in equal length...")
                self.has_metadata = True
                self.has_neighbors = True
        else:
            logger.debug("KDB instance has no metadata or neighbor information...")
            self.kmer_utils = kmer_utils.Utils(self.header['k'])
        if make_neighbors is True and self.header["metadata"] is True:
            if self.has_neighbors is False:
                logger.debug("kdb.fileutil.KDB.__init__() is creating neighbor information")
            else:
                logger.warning("kdb.fileutil.KDB.__init__() is re-creating neighbor information")                
            self._make_neighbors()
        elif self.has_neighbors is True:
            logger.debug("kdb.fileutil.KDB.__init__() was initialized with a neighbor array")


    def _make_neighbors(self):
        if self.has_neighbors is True:
            logger.warn("kdb.fileutil.KDB._make_neighbors called on top of an existing neighbors array")
            return
        else:
            logger.debug("Generating the neighbor array...")
            for x in range(4**self.header['k']):
                self.neighbors.append(None)
            for idx, c in enumerate(self.profile):
                if c > 0:
                    # neighbors = list(filter(lambda x: self.profile[x] > 0, kmer_utils.get_neighbors(idx, self.header['k'])))
                    # for y in neighbors:
                    #     logger.debug("\t".join(list(map(str, [idx, c, y, kmer_utils.binary_to_sequence(y, self.header['k']), self.profile[y]]))))
                    
                    self.neighbors[idx] = [
                        {
                            "id": y,
                            "count_delta": abs(self.profile[y] - c) # The difference between k-mer counts
                            # A markov probability delta (0-1) is a one-to-many between neighbors, node centric
                            #@"markov": vector_ops.markov_probability(profile)
                        } for y in list(filter(lambda x: self.profile[x] > 0, self.kmer_utils.get_neighbors(idx)))]
            self.has_neighbors = True

    def _make_metadata(self):
        if self.has_metadata is True:
            logger.warn("kdb.fileutil.KDB._make_metadata called on top of an existing neighbors array")
            return
        else:
            logger.debug("Generating the metadata array...")
            for x in range(4**self.header['k']):
                self.metadata.append(None)
            for idx, c in enumerate(self.profile):
                if c > 0:
                    metadata = {
                        "tags": [],
                        "bools": [],
                        "floats": [],
                        "ints": [],
                        "neighbors": [
                            {
                                "id": x,
                                "count_delta": abs(self.profile[idx] - c) # The difference between k-mer counts
                                # A markov probability delta (0-1) is a one-to-many between neighbors, node centric
                                #@"markov": vector_ops.markov_probability(profile)
                            } for x in list(filter(lambda x: self.profile[x] > 0, self.kmer_utils.get_neighbors(idx)))]
                    }
                    if self.has_neighbors is True:
                        metadata["neighbors"] = self.neighbors[idx]
                    self.metadata[idx] = metadata
            self.has_metadata = True

            
    # def _regenerate_header(self):
    #     self.total_kmers = functools.reduce(lambda a,b: a+b, self.profile)
    #     self.nullomers = self.profile.count(0)
    #     self.unique_kmers = 4**self.k - self.nullomers
    #     self.total_reads = i
    #     self.header["total_kmers"]  = self.total_kmers
    #     self.header["nullomers"]    = self.nullomers
    #     self.header["unique_kmers"] = self.unique_kmers
    #     self.total_reads
                    
    def add_profile(self, x: array.array):
        if type(x) is not array.array:
            raise TypeError("kdb.fileutil.KDB.add_profile expects an array as its argument")
        elif 4**self.header['k'] != len(x):
            raise TypeError("kdb.fileutil.KDB.add_profile expects the arrays to have equal length")
        for i in range(len(x)):
            self.profile[i] += x[i]
        #self._regenerate_header()

        
            
                    
class KDBReader:
    def __init__(self, fileobj: io.IOBase, include_neighbors=False, include_metadata=False):
        if not isinstance(fileobj, io.IOBase):
            raise TypeError("kdb.fileutil.KDBReader expects a file object as its first positional argument")
        elif type(include_neighbors) is not bool:
            raise TypeError("kdb.fileutil.KDBReader expects the keyword argument 'include_neighbors' to be a bool")
        elif type(include_neighbors) is not bool:
            raise TypeError("kdb.fileutil.KDBReader expects the keyword argument 'include_metadata' to be a bool")
        self.include_neighbors = include_neighbors
        self.include_metadata  = include_metadata
        self.has_metadata = False
        self.loaded       = False
        self.filehandle   = fileobj
        self.filepath     = self.filehandle.name
        self.profile      = array.array('L')#array.array('H')
        self.metadata     = []
        self.neighbors    = []
        basename, ext     = os.path.splitext(self.filepath)
        # Extension validation
        if ext == '':
            raise IOError("kdb.fileutil.KDBReader could not determine the extention of '{0}'".format(self.filepath))
        elif ext != '.bgzf' and ext != ".kdb":
            raise IOError("kdb.file.get_header can only parse .bgzf format '.kdb' files.")
        # Offset handling
        self._bgzfrdr = bgzf.BgzfReader(fileobj=self.filehandle)
        self._offsets = []
        # # # # # # # # # # # # #
        # Parse header
        # # # # # # # # # # # # #

        # Ensure the pointers are at the start of the file
        self._bgzfrdr.seek(0) 
        self.filehandle.seek(0)
        for values in bgzf.BgzfBlocks(self.filehandle):
            #logger.debug("Raw start %i, raw length %i, data start %i, data length %i" % values)
            self._offsets.append(values) # raw start, raw length, data start, data length
        if len(self._offsets) == 0:
            raise IOError("kdb.fileutil.KDBReader opened an empty file")
        # Ensure the pointers are at the start of the file
        self._bgzfrdr.seek(0) 
        self.filehandle.seek(0)
        # Read data-length bytes from the zeroth(header) block
        header=self._bgzfrdr.read(self._offsets[0][3]) 
        logger.info("Read KDB header:\n{0}".format(header))
        data=yaml.safe_load(header)
        if type(data) is not dict:
            raise TypeError("kdb.fileutil.KDBReader the header of '{}' was not parseable YAML formatted".format(self.filepath))
        else:
            try:
                jsonschema.validate(instance=data, schema=config.header_schema)
                data["body_start"] = self._offsets[0][3]+1 # Ensure header/body offset information is stored in the header while in memory
                self.header = data
                self.k = self.header['k']
            except jsonschema.ValidationError as e:
                logger.debug(e)
                logger.error("kdb.fileutil.KDBReader couldn't validate the header YAML")
                raise e
            if self.include_neighbors:
                if self.include_metadata:
                    for x in range(4**self.k):
                        self.profile.append(0)
                        self.neighbors.append(None)
                        self.metadata.append(None)
                else:
                    for x in range(4**self.k):
                        self.profile.append(0)
                        self.neighbors.append(None)
            elif self.include_metadata:
                for x in range(4**self.k):
                    self.profile.append(0)
                    self.metadata.append(None)
            else:
                for x in range(4**self.k):
                    self.profile.append(0)




            
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()
        
    def read_profile(self):
        """ Read profile from the .kdb file. Stores attribute in the .profile attribute.
        Does not mutate the KDBReader object, just reads data from the file. If the file header doesn't validate, it shouldn't get to this point.
        
        :returns: 
        :rtype: NoneType

        """
        logger.info("Reading profile from '{}'...".format(self.filepath))
        offsets = copy.deepcopy(self._offsets)
        # Omit the header
        offsets.pop(0)
        # Ensure the pointers are at the start of the file
        self._bgzfrdr.seek(0) 
        #self.filehandle.seek(0)

        # Read data-length bytes from the zeroth(header) block
        logger.debug("Reading blocks from '{}'...".format(self.filepath))
        for raw_start, raw_length, data_start, data_length in offsets:
            if data_length > 0:
                #logger.debug("Raw start: {0}, Raw length: {1}, Data start: {2}, Data length: {3}".format(raw_start, raw_length, data_start, data_length))
                offset = bgzf.make_virtual_offset(raw_start, 0)
                self._bgzfrdr.seek(offset)
                #self.filehandle.seek(data_start)
                data = self._bgzfrdr.read(data_length)

                for line in data.split("\n"):
                    if line != '':
                        # FIXME
                        # This part should include more data definition.
                        # Type validation, neighbor validation, parameterized on
                        # data types, data mutability constraints, eventually subject to the index.

                        
                        ################
                        # DATA PARSING
                        ################
                        # JSON parsing
                        idx, count, metadata = line.split("\t")
                        idx, count = int(idx), int(count)
                        metadata = yaml.safe_load(metadata) # Exceptions shouldn't be caught
                        if len(metadata) != 0:
                            try:
                                jsonschema.validate(instance=metadata, schema=config.metadata_schema)

                                self.has_metadata = True
                                if "neighbors" in metadata:
                                    neighbors = metadata["neighbors"]
                                else:
                                    logger.warning("kdb.fileutil.KDBReader couldn't validate a k-mer's metadata YAML : no key 'neighbor'")
                                    neighbors = None
                                    self.has_neighbors = False
                            except jsonschema.ValidationError as e:
                                logger.debug(e)
                                logger.error("kdb.fileutil.KDBReader couldn't validate a k-mer's metadata YAML")
                                raise e
                        else:
                            metadata = None
                            neighbors = None
                        
                        # Single-line parsing
                        #idx, count, neighbors = line.split("\t")
                        #self.profile[int(idx)] = int(count)

                        self.profile[idx] = count
                        if self.include_metadata:
                            self.metadata[idx] = metadata
                        if self.include_neighbors:
                            self.neighbors[idx] = neighbors
                self.loaded = True
                    
        
class KDBWriter:
    def __init__(self, fileobj: io.IOBase, kdb: KDB):
        if not isinstance(fileobj, io.IOBase):
            raise TypeError("kdb.fileutil.KDBWriter expects a file object as its first positional argument")
        elif not isinstance(kdb, KDB):
            raise TypeError("kdb.fileutil.KDBWriter expects a KDB objects as its second positional argument")
        self.kdb = kdb
        self.filehandle = bgzf.BgzfWriter(fileobj=fileobj, compresslevel=9)
        self.filepath = fileobj.name
        
        # Write header
        header = kdb.header
        # Remove header/body offset information
        if 'body_start' in header:
            header.pop('body_start')
        self.filehandle._write_block(bgzf._as_bytes(yaml.dump(header)))
        self.filehandle.flush()

            
    def __exit__(self, exc_type, exc_value, traceback):
        self.filehandle.close()
        
    def write_profile(self, include_metadata):
        """ Write the profile to the file defined during KDBWriter __init__. This is just the wrapper for the IO operations.
        :param include_metadata: Whether or not to write the k-mer metadata to the file
        :type include_metadata: bool
        :returns: 
        :rtype: NoneType
        
        """
        if type(include_metadata) is not bool:
            raise TypeError("kdb.fileutil.KDBWriter.write_profile expects a bool as its only positional argument")
            
        logger.info("Writing profile to  '{}'...".format(self.filepath))
        # Write body
        i=0

        text=''
        size = 0
        for idx, c in enumerate(self.kdb.profile):
            #logger.debug("Writing k-mer to file: {0} - {1}".format(idx, c))
            node = {
                "id": idx,
                "count": c,
                "metadata": {
                    "tags": [],
                    "bools": [],
                    "floats": [],
                    "ints": []
                }
            }

            if c > 0:
                if include_metadata and self.kdb.has_metadata: # Use existing metadata
                    node["metadata"] = self.kdb.metadata[idx]
                elif include_metadata and self.kdb.has_neighbors: # Use existing neighbor
                    node["metadata"]["neighbors"] = self.kdb.neighbors[idx]
                elif include_metadata: # Create the k-mer neighbor metadata from scratch
                    neighbors = list(filter(lambda x: self.kdb.profile[x] > 0, self.kmer_utils.get_neighbors(idx)))
                    node["metadata"]["neighbors"] = [
                        {
                            "id": x,
                            "count_delta": abs(self.kdb.profile[idx] - c) # The difference between k-mer counts
                            # A markov probability delta (0-1) is a one-to-many between neighbors, node centric
                            #@"markov": vector_ops.markov_probability(profile)
                        } for x in neighbors
                    ]
            else:
                node["metadata"] = {}
            newline = "{0}\t{1}\t{2}\n".format(node["id"], node["count"], json.dumps(node["metadata"]))
            #newline = json.dumps(node)
            size += bgzf._as_bytes(newline).__sizeof__()
            i+=1
            if i%1000 == 0:
                sys.stderr.write("\r")
                sys.stderr.write("Wrote {} lines to the file...".format(i))
            #logger.debug("Line:\n{0}Size in bytes: {1}".format(newline, size))
            if size > 65534:
                self.filehandle._write_block(bgzf._as_bytes(text)) # should this be _write_block?
                self.filehandle._buffer = b""
                self.filehandle._handle.flush()
                #self.filehandle.flush()
                text = ''
                size = 0
            text += newline

        if text != '': # Final flush of text buffer
            sys.stderr.write("\r")
            sys.stderr.write("Wrote {} lines to the file...\n".format(i))
            logger.info("Writing final block...")
            self.filehandle._write_block(bgzf._as_bytes(text))
            self.filehandle._buffer = b""
            self.filehandle._handle.flush()

            #self.filehandle.flush()

        if i != self.kdb.num_kmers:
            raise IOError("Wrote {0} lines (of {1} k-mers) to the file...".format(i, self.kdb.num_kmers))
