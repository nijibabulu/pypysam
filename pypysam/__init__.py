#! /usr/bin/env python

import os
import types
import string
import itertools
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
        self.cur_buffer_size = len(self.chunk)
        self.buffer_pos = 0
    def tell(self):
        return self.pos

    def read(self):
        self.pos += 1
        if self.buffer_pos == self.cur_buffer_size:
            self._next_chunk()
            if not self.chunk:
                return ''
        self.buffer_pos += 1
        return self.chunk[self.buffer_pos-1]

    def seek(self,pos):
        os.lseek(self.fileno,pos,os.SEEK_SET)
        self.buffer_pos = self.cur_buffer_size = 0

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
            if c.isalpha():
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



class AlignmentReference(object):
    def __init__(self,tid,name,length):
        self.tid = tid
        self.name = name
        self.length = length
    def __cmp__(self,other):
        return cmp(self.tid,other.tid) or cmp(self.name,other.name) or cmp(self.length,other.length)

class AlignmentFileBase(object):
    def __init__(self,filepath,mode='r',template=None,header=None,**kwargs):
        self.filepath = filepath
        self.mode = mode
        self.is_bam = False
        self.is_cram = False
        self.is_sam = False
        self.write_mode = False
        self.reference_table = collections.OrderedDict()
        self.reference_names = []
        if mode[0] == 'r':
            self.read_header()
        elif mode[0] == 'w':
            if template is not None:
                self.header_from_template(template)
            elif header is not None:
                self.header_from_object(header)
            self.write_header()
            self.write_mode = True
        else:
            raise ValueError, 'Unknown mode %s' % mode

class AlignedSegment(object):
    _cigar_ops = ['M','I','D','N','S','H','P','=','X']
    def __init__(self,aln_file,tid,pos,mapq,flag,mate_tid,mate_pos,
                 tlen,read_name,cigar,seq,qual,tags,block=None):
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
        self.tags = tags
        self.block = block
        self._block_format = aln_file._magic
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
        self._value = value

    @property
    def value(self):
        return self._value

class SamTag(Tag):
    _tag_sam_types = 'cCsSiIAf'
    _tag_pystruct_types = 'bBhHiIcf'
    _tag_trans = string.maketrans(_tag_sam_types,_tag_pystruct_types)

    def __init__(self,tag,val_type,value=None,value_block=None):
        super(SamTag,self).__init__(tag,val_type,value)
        self._value_block = value_block
        if value is None and value_block is None:
            raise ValueError, 'Tag requires either a value or a value block'

    @property
    def value(self):
        if self._value is None:
            if self.val_type == 'Z':
                self._value = self._value_block
            elif self.val_type in self._tag_sam_types:
                self._value = struct.unpack(
                    self.val_type.translate(self._tag_trans),self._value_block)[0]
            else:
                raise NotImplementedError, 'Tag type %s not implemented' % val_type
        return self._value

    @property
    def value_block(self):
        if self._value_block is None:
            self._value_block = struct.pack(
                tag.val_type.translate(self._tag_trans),self._value)
        return self._value_block


class SamFile(AlignmentFileBase):
    _magic = b"\x40\x48\x44\x20"
    def __init__(self):
        raise NotImplementedError, 'SAM file handling is not implemented yet'

           
class CramFile(AlignmentFileBase):
    _magic = b"\x43\x52\x41\x4d"
    def __init__(self):
        raise NotImplementedError, 'Cram file handling is not implemented yet'

def AlignmentFile(filepath,mode='r',**kwargs):
    filepath = filepath
    fh = open(filepath,mode)
    if mode[0] == 'r':
        magic = fh.read(4)
        magic_map = {cls._magic:cls for cls in [BamFile,CramFile,SamFile]}
        if magic not in magic_map:
            raise ValueError, 'Unknown format for file %s' % filepath
        return magic_map[magic](filepath,mode,**kwargs)
    elif mode[0] == 'w':
        if 'format' not in kwargs:
            raise ValueError, 'Specify an output format for writing'
        fmt = kwargs['format']
        format_map = {cls.__name__[:len(fmt)].lower():cls 
                      for cls in [BamFile,CramFile,SamFile]}
        if fmt not in format_map:
            raise ValueError, 'Unknown format %s' % fmt
        return format_map[fmt.lower()](filepath,mode,**kwargs)


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

# calculate bin given an alignment covering [beg,end) (zero-based, half-closed-half-open) 
# from SAM format specification
def reg2bin(beg, end):
    end -= 1
    if beg>>14 == end>>14: return ((1<<15)-1)/7 + (beg>>14);
    if beg>>17 == end>>17: return ((1<<12)-1)/7 + (beg>>17);
    if beg>>20 == end>>20: return ((1<<9)-1)/7 + (beg>>20);
    if beg>>23 == end>>23: return ((1<<6)-1)/7 + (beg>>23);
    if beg>>26 == end>>26: return ((1<<3)-1)/7 + (beg>>26);
    return 0;

# calculate the list of bins that may overlap with region [beg,end) (zero-based) 
# from SAM format specification
MAX_BIN=(((1<<18)-1)/7)
def reg2bins(beg, end):
    end -= 1
    bin_list = [0]

    for k in range(1 + (beg>>26), 2 + (end>>26)): bin_list.append(k)
    for k in range(9 + (beg>>23), 10 + (end>>23)): bin_list.append(k)
    for k in range(73 + (beg>>20), 74 + (end>>20)): bin_list.append(k)
    for k in range(585 + (beg>>17), 586 + (end>>17)): bin_list.append(k)
    for k in range(4681+ (beg>>14), 4682 + (end>>14)): bin_list.append(k)

    return bin_list

def _bam_cigar_type(val):
    return _bam_cigar_type.c>>((val)<<1)&3
_bam_cigar_type.c = 0x3C1A7

class BamFile(AlignmentFileBase):
    _seq_dec = '=ACMGRSVTWYHKDBN' 
    _seq_enc = {x:i for i,x in enumerate(list(_seq_dec))}
    _magic = bgzf._bgzf_magic
    _header_magic = 'BAM\1'

    _tag_sam_types = 'cCsSiIAf'
    _tag_pystruct_types = 'bBhHiIcf'
    _tag_trans = string.maketrans(_tag_sam_types,_tag_pystruct_types)
    _tag_sizes = {'c':1, 'C':1, 's':2, 'S':2, 'i': 4, 'I': 4, 'A': 1, 'f': 4}

    def __iter__(self):
        return self

    def write(self,aln):
        if not self.write_mode:
            raise ValueError, 'Not opened for writing'
        if aln.block is not None and aln._block_format == self._magic:
            block = aln.block
        else:
            seq_len = len(aln.seq)
            block = struct.pack(
                'iiiiiiii',
                aln.tid,
                aln.pos,
                reg2bin(aln.pos,aln.ref_alignment_end) << 16 | 
                    aln.mapq << 8 | (len(aln.qname)+1),
                aln.flag << 16 | len(aln.cigar),
                seq_len,
                aln.mate_tid,
                aln.mate_pos,
                aln.tlen) + aln.qname + '\x00'

            for op,op_len in aln.cigar:
                block += struct.pack('i',op_len << 4 | op)
            for c1,c2 in itertools.izip(itertools.islice(aln.seq,0,None,2),
                                        itertools.islice(aln.seq,1,None,2)):
                block += struct.pack('B', self._seq_enc[c1]<<4 | self._seq_enc[c2])
            if seq_len % 2:
                block += struct.pack('B',self._seq_enc[aln.seq[-1]] << 4)
            #block += ''.join(chr(c) for c in aln.qual)
            block += struct.pack('b'*seq_len,*aln.qual)
            for tag in aln.tags:
                block += tag.tag+tag.val_type
                if tag.val_type in self._tag_sam_types:
                    block += tag.value_block
                elif tag.val_type == 'Z':
                    block += tag.value + '\0'
                else:
                    raise ValueError, 'Tag type %s not supported' % tag.val_type

        self.fh.write(struct.pack('i',len(block)) + block)

    def next(self):
        block_size_block = self.fh.read(4)
        if len(block_size_block) == 0:
            raise StopIteration
        block_size = struct.unpack('i',block_size_block)[0]
        assert block_size > 0

        block = self.fh.read(block_size)

        refID,pos,bin_mq_nl,flag_nc,l_seq,next_refID,next_pos,tlen =\
            struct.unpack('iiIIiiii',block[:32])

        mapq = (bin_mq_nl >> 8) & 0xff
        flag = flag_nc >> 16
        read_name = block[32:].split('\x00')[0]
        block = block[33+len(read_name):]

        n_cigar_op = flag_nc & 0xffff
        cigar_end = n_cigar_op*4
        raw_cigar = struct.unpack('I'*n_cigar_op,block[:cigar_end])
        cigar = tuple((x&0xf,x>>4) for x in raw_cigar)

        p = cigar_end
        seq_end = p+(l_seq+1)/2
        enc_seq = struct.unpack('B'*((l_seq+1)/2),block[p:seq_end])
        seq = ''
        for c in enc_seq:
            seq += self._seq_dec[c>>4]
            if len(seq) < l_seq:
                seq += self._seq_dec[c & 0xf]

        p = seq_end
        qual_end = p+l_seq

        qual = struct.unpack('b'*l_seq,block[p:qual_end])

        remainder = block_size - 32 - len(read_name) - 1 - 4*n_cigar_op -\
                (l_seq+1)/2 - l_seq
        tag_block = block[p+l_seq:]

        i = 0
        tags = []
        while i < len(tag_block)-1:
            tag = tag_block[i:i+2]
            val_type = tag_block[i+2]
            if val_type in self._tag_sam_types:
                value_block = tag_block[i+3:i+3+self._tag_sizes[val_type]]
            elif val_type == 'Z':
                value_block = tag_block[i+3:].split('\x00')[0]
            else:
                raise ValueError, 'Tag type %s not supported' % tag.val_type
            i += 3+len(value_block)
            tags.append(SamTag(tag,val_type,value_block=value_block))

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
        self.header_text = self.fh.read(l_text)
        n_ref = struct.unpack('i',self.fh.read(4))[0]
        for tid in range(n_ref):
            l_name = struct.unpack('i',self.fh.read(4))[0]
            name = self.fh.read(l_name)
            name = name[:-1]
            l_ref = struct.unpack('i',self.fh.read(4))[0]
            self.reference_table[name] = AlignmentReference(tid,name,l_ref)
            self.reference_names.append(name)

    def header_from_template(self,template):
        self.header_text = template.header_text
        self.reference_table = template.reference_table

    # XXX: make a single object for the header instead of text and table
    def header_from_object(self,object):
        raise NotImplementedError, 'Headers from objects not supported yet'

    def write_header(self):
        if not hasattr(self,'header_text'):
            raise ValueError, 'Cannot write a header when none is defined'
        self.fh = bgzf.BgzfWriter(self.filepath,self.mode,compresslevel=1)
        block = self._header_magic
        block += struct.pack('i',len(self.header_text)+1)
        block += self.header_text + '\0'
        block += struct.pack('i',len(self.reference_table))
        #print len(self.reference_table)
        for name,aref in self.reference_table.items():
            block += struct.pack('i',len(name)+1)
            block += name + '\0'
            block += struct.pack('i',aref.length)

        self.fh.write(block)
            #print len(name)+1,name,aref.length

 


def main():
    import sys

    #if len(sys.argv) > 1:
        #fa = FastaFile(sys.argv[1])
        #raise SystemExit
    if len(sys.argv) > 1:
        af = AlignmentFile(sys.argv[1])
        #for aln in af:
            #print aln.seq,aln.cigar,aln.cigarstring
        if len(sys.argv) > 2:
            out = AlignmentFile(sys.argv[2],'wb',format='bam',template=af)
            for aln in af:
                out.write(aln)
            out.fh.close()
            #af2 = AlignmentFile(sys.argv[2])

            #for (k1,v1),(k2,v2) in zip(af.reference_table.items(),
                                       #af2.reference_table.items()):
                #print k1,k2,k1==k2,v1==v2
            #pprint.pprint(af.reference_table)
            #pprint.pprint(af2.reference_table)



if __name__ == '__main__':
    import cProfile
    import pstats
    cProfile.run('main()','tmpstats')
    p = pstats.Stats('tmpstats')
    p.strip_dirs().sort_stats('cumtime').print_callees()
    #main()
