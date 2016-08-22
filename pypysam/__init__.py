#! /usr/bin/env python

import os
import types
import pprint
import struct
import collections
import bgzf

class BufferedOSRead(object):
    ''' using os.read in pypy is faster than BufferedIO or default IO'''
    def __init__(self,path,buffer_size=16384):
        self.path = path
        self.f = open(self.path)
        self.fileno = self.f.fileno()
        self.buffer_size = buffer_size
        self.pos = 0
        self._next_chunk()
    def _next_chunk(self):
        self.chunk = os.read(self.fileno,self.buffer_size)
        self.buffer_pos = 0
    def tell(self):
        return self.pos

    def read(self):
        self.pos += 1
        if self.buffer_pos == self.buffer_size:
            self._next_chunk()
            if not self.chunk:
                return None
        self.buffer_pos += 1
        return self.chunk[self.buffer_pos-1]

    def seek(self,pos):
        os.lseek(self.fileno,pos,os.SEEK_SET)

    def readline(self):
        s = ''
        while True:
            s += self.read()
            if s[-1] == '\n':
                return s

def parse_region(region_str):
    try:
        colon_sep = region_str.split(':')
        if len(colon_sep) == 2:
            start,end = [int(x) for x in colon_sep[1].split('-')]
        elif len(colon_sep) == 3:
            start,end = [int(x) for x in colon_sep[1:]]
        return colon_sep[0],start,end
    except: 
        raise ValueError, 'Improperly fomatted region %s' % region_str

class IndexedSequence(object):
    __slots__ = ('name','length','position','line_len','pad_line_len')
    def __init__(self,name,length,position,line_len,pad_line_len):
        self.name = name
        self.length = int(length)
        self.position = int(position)
        self.line_len = int(line_len)
        self.pad_line_len = int(pad_line_len)

class FastaFile(object):
    def __init__(self,path):
        self.path = path
        self.idx_path = self.path + '.fai'
        self.f = BufferedOSRead(path)
        if not os.path.exists(self.idx_path):
            self._build_index()
        self._read_index()

    def fetch(self,reference=None,start=None,end=None,region=None):
        if region is not None:
            reference,start,end = parse_region(region)
        elif reference is None:
            raise ValueError, 'fetch() requires either a reference or a region'
        ind_seq = self.seqs[reference]
        fstart = ind_seq.position
        if start is not None:
            fstart += start / ind_seq.line_len * ind_seq.pad_line_len + \
                start % ind_seq.line_len
        self.f.seek(fstart)
        ret_seq = ''
        while True:
            c = self.f.read()
            if c == '>' or c == '':
                break
            if not c.isspace():
                ret_seq += c
            if end is not None and len(ret_seq) >= end-start:
                break
        return ret_seq

    def _read_index(self):
        self.seqs = {}
        for line in open(self.idx_path):
            fields = line.split()
            self.seqs[fields[0]] = IndexedSequence(*fields)

    def _build_index(self):
        self.f = open(self.path)
        self.fileno = self.f.fileno()
        self.br = BufferedOSRead(self.path)
        out = open(self.idx_path,'w')

        class _FaState():
            def __init__(self):
                self.name = self.tell = self.line_len = self.pad = -1
                self.cur_line_len = self.cur_pad = self.length = 0
                self.at_end = False
        s = _FaState()

        counter = 0
        line_no = 0
        while True:
            c = self.br.read()
            #c = self.f.read(1)

            if c == '\n': 
                line_no += 1
            if not c: # error handling not right here!
                out.write('\t'.join([str(x) for x in 
                                     [s.name,s.length,s.tell,
                                      s.line_len,s.line_len+s.pad]])+'\n')
                out.close()
                break
            if c == '>':
                if s.name != -1:
                    out.write('\t'.join([str(x) for x in 
                                         [s.name,s.length,s.tell,
                                          s.line_len,s.line_len+s.pad]])+'\n')
                s = _FaState()
                #s.name = self.f.readline().split()[0]
                s.name = self.br.readline().split()[0]
                line_no += 1
                s.tell = self.br.tell()
            elif s.at_end:
                raise IOError, 'Fasta file must have same line lengths'
            elif c.isspace():
                if s.cur_pad == 0:
                    s.cur_pad = 1
                    if s.line_len == -1:
                        s.line_len = s.cur_line_len
                    elif s.line_len != s.cur_line_len:
                        # odd line length must indicate end of sequence
                        s.at_end = True 
                    s.cur_line_len = 0
                else:
                    s.cur_pad += 1
            else:
                if s.cur_pad > 0:
                    if s.pad == -1:
                        s.pad = s.cur_pad
                    elif s.cur_pad != s.pad:
                        raise IOError, \
                                'Odd number of line feeds in line %d' % line_no
                    s.cur_pad = 0
                s.cur_line_len += 1
                s.length += 1


_cram_magic = b"\x43\x52\x41\x4d"
_sam_magic = b"\x40\x48\x44\x20"

class AlignmentReference(object):
    def __init__(self,tid,name,length):
        self.tid = tid
        self.name = name
        self.length = length

class AlignmentFileBase(object):
    def __init__(self,filepath):
        self.filepath = filepath
        self.is_bam = False
        self.is_cram = False
        self.is_sam = False
        self.reference_table = collections.OrderedDict()
        self.reference_names = []
        self.read_header()

class cached_property(object):
    """ A property that is only computed once per instance and then replaces
        itself with an ordinary attribute. Deleting the attribute resets the
        property.

        Source: https://github.com/bottlepy/bottle/commit/fa7733e075da0d790d809aa3d2f53071897e6f76
        """

    def __init__(self, func):
        self.__doc__ = getattr(func, '__doc__')
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return self
        value = obj.__dict__[self.func.__name__] = self.func(obj)
        return value

class AlignedSegment(object):
    _cigar_ops = ['M','I','D','N','S','H','P','=','X']
    def __init__(self,aln_file,tid,pos,mapq,flag,mate_tid,mate_pos,
                 tlen,read_name,cigar,seq,qual,tags):
        self.qname = read_name
        self.tid = tid
        self.pos = pos
        self.mapq = mapq
        self.flag = flag
        self.seq = seq
        self.mate_tid = mate_tid
        self.mate_pos = mate_pos
        self.tlen = tlen
        self.read_name = read_name
        self.cigar = cigar
        self.qual = qual
        self.tags = {t.tag:t for t in tags}
        self._aln_file = aln_file

    def has_tag(self,tag_name):
        return tag_name in self.tags
    @property
    def query_length(self):
        return len(self.seq)

    @property
    def reference_name(self):
        return self._aln_file.reference_names[self.tid]

    @property
    def cigarstring(self):
        return str(
            ''.join(''.join((str(c[1]), self.__class__._cigar_ops[c[0]])) \
                       for c in self.cigar))

    @property
    def query_sequence(self):
        return self.seq

    @property
    def query_alignment_start(self):
        start = 0
        if self.cigar[0][0] == 4:
            start = self.cigar[0][1]
        return start


    @property
    def query_alignment_end(self):
        end = 0
        for op,length in self.cigar:
            if _bam_cigar_type(op) & 1:
                end += length
        return end

    @property
    def ref_alignment_end(self):
        end = self.pos
        for op,length in self.cigar:
            if _bam_cigar_type(op) & 2:
                end += length
        return end

    @property
    def aend(self):
        return self.ref_alignment_end

    @property
    def query_alignment_sequence(self):
        return self.query_sequence[self.query_alignment_start:
                                   self.query_alignment_end]

    @property
    def is_proper_pair(self):
        return self.flag & BamFlag.PROPER_PAIR

    @property
    def is_duplicate(self):
        return self.flag & BamFlag.DUP

    @property
    def is_paired(self):
        return self.flag & BamFlag.PAIRED

    @property
    def is_qcfail(self):
        return self.flag & BamFlag.QCFAIL

    @property
    def is_read1(self):
        return self.flag & BamFlag.IS_READ1

    @property
    def is_read2(self):
        return self.flag & BamFlag.IS_READ2

    @property
    def is_reverse(self):
        return self.flag & BamFlag.REVERSE

    @property 
    def is_secondary(self):
        return self.flag & BamFlag.SECONDARY

    @property
    def is_supplementary(self):
        return self.flag & BamFlag.SUPPLEMENTARY

    @property
    def is_unmapped(self):
        return self.flag & BamFlag.UNMAP

class Tag(object):
    def __init__(self,tag,val_type,value):
        self.tag = tag
        self.val_type = val_type
        self.value = value

class SamFile(AlignmentFileBase):
    def read_header(cls,obj):
        raise NotImplementedError, 'Sam parsing is not implemented yet'

           
class CramFile(AlignmentFileBase):
    def read_header(cls,obj):
        raise NotImplementedError, 'Cram parsing is not implemented yet'
    pass

def AlignmentFile(filepath):
    filepath = filepath
    fh = open(filepath)
    magic = fh.read(4)
    if magic == bgzf._bgzf_magic:
        return BamFile(filepath)
    elif magic == _cram_magic:
        return CramFile(filepath)
    elif magic == _sam_magic:
        return SamFile(filepath)
    else:
        raise ValueError, 'Unknown format for file %s' % filepath


class BamFlag(object):
    PAIRED = 0x1
    PROPER_PAIR = 0x2
    UNMAP = 0x4
    MUNMAP = 0x8
    REVERSE = 0x10
    MREVERSE = 0x20
    READ1 = 0x40
    READ2 = 0x80
    SECONDARY = 0x100
    QCFAIL = 0x200
    DUP = 0x400
    SUPPLEMENTARY = 0x800

def _bam_cigar_type(val):
    return _bam_cigar_type.c>>((val)<<1)&3
_bam_cigar_type.c = 0x3C1A7

class BamFile(AlignmentFileBase):
    _seq_enc = '=ACMGRSVTWYHKDBN' 

    def __iter__(self):
        return self

    def next(self):
        block_size_block = self.fh.read(4)
        if len(block_size_block) == 0:
            raise StopIteration
        block_size = struct.unpack('i',block_size_block)[0]
        assert block_size > 0

        block = self.fh.read(block_size)

        refID = struct.unpack('i',block[0:4])[0]
        pos = struct.unpack('i',block[4:8])[0]
        # note don't know about bin yet--probably need to store or compute this?
        bin_mq_nl = struct.unpack('I',block[8:12])[0]
        mapq = (bin_mq_nl >> 8) & 0xff
        flag_nc = struct.unpack('i',block[12:16])[0]
        flag = flag_nc >> 16
        l_seq = struct.unpack('i',block[16:20])[0]
        next_refID = struct.unpack('i',block[20:24])[0]
        next_pos = struct.unpack('i',block[24:28])[0]
        tlen = struct.unpack('i',block[28:32])[0]
        read_name = block[32:].split('\x00')[0]
        block = block[33+len(read_name):]
        #while True:
            #c = self.fh.read(1)
            #if ord(c) == 0:
                #break
            #read_name += c

        n_cigar_op = flag_nc & 0xffff
        raw_cigar = [struct.unpack('I',block[i*4:i*4+4])[0] 
                 for i in range(n_cigar_op)]
        cigar = tuple((x&0xf,x>>4) for x in raw_cigar)

        p = n_cigar_op*4
        enc_seq = [struct.unpack('B',block[p + i])[0] 
                   for i in range((l_seq+1)/2)]
        seq = ''
        for c in enc_seq:
            seq += self._seq_enc[c>>4]
            if len(seq) < l_seq:
                seq += self._seq_enc[c & 0xf]

        p += (l_seq+1)/2

        qual = [struct.unpack('b',block[p+i])[0]
                for i in range(l_seq)]

        remainder = block_size - 32 - len(read_name) - 1 - 4*n_cigar_op -\
                (l_seq+1)/2 - l_seq
        tag_block = block[p+l_seq:]
        #tag_block = self.fh.read(remainder)
        i = 0
        tags = []
        while i < len(tag_block)-1:
            tag = tag_block[i:i+2]
            val_type = tag_block[i+2]
            if val_type == 'Z':
                i += 3
                value = ''
                while True:
                    c = tag_block[i]
                    i += 1
                    if ord(c) == 0:
                        break
                    value += c
            elif val_type == 'c' or val_type == 'C':
                value = struct.unpack(chr(ord(val_type)-1),tag_block[i+3])[0]
                i += 4
            elif val_type == 's' or val_type == 'S':
                value = struct.unpack(chr(ord(val_type)-11),tag_block[i+3:i+5])[0]
                i += 5
            elif val_type == 'i' or val_type == 'I':
                value = struct.unpack(val_type,tag_block[i+3:i+7])[0]
                i += 7
            elif val_type == 'A':
                value = tag_block[i+3]
                i += 4
            elif val_type == 'f':
                value = struct.unpack(val_type,tag_block[i+3:i+7])[0]
                i += 7
            else:
                raise NotImplementedError, 'Tag type %s not implemented' % val_type
            tags.append(Tag(tag,val_type,value))
            return AlignedSegment(self,refID,pos,mapq,flag,next_refID,next_pos,
                                  tlen,read_name,cigar,seq,qual,tags)
            
        '''
        Z - string
        A - character
        i - signed integer
        f - single float
        H - byte array in hex format
        B - [cCsSiIf] integer or numeric array
            bytes 5-8 (0xf) store the number of elements in the array

        note cCsSiI are all "i" in SAM format. have to convert this back in bam
        format based on the size of the integer (i guess)
        '''
        

    def read_header(self):
        self.fh = bgzf.open(self.filepath)
        magic = struct.unpack('i',self.fh.read(4))[0]
        l_text = struct.unpack('i',self.fh.read(4))[0]
        text = self.fh.read(l_text)
        n_ref = struct.unpack('i',self.fh.read(4))[0]
        for tid in range(n_ref):
            l_name = struct.unpack('i',self.fh.read(4))[0]
            name = self.fh.read(l_name)
            name = name[:-1]
            l_ref = struct.unpack('i',self.fh.read(4))[0]
            self.reference_table[name] = AlignmentReference(tid,name,l_ref)
            self.reference_names.append(name)
 


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        fa = FastaFile(sys.argv[1])
        raise SystemExit
    if len(sys.argv) > 1:
        af = AlignmentFile(sys.argv[1])
        for aln in af:
            print aln.seq,aln.cigar,aln.cigarstring



