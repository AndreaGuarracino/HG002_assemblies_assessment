#!/usr/bin/env python
## category General
## desc Postprocesses a BAM file to remove all clipping from reads and alignments
## experimental
'''
Postprocesses a BAM file to remove all clipping from reads and alignments.

Hard clipping is removed from the alignment. Soft clipping is removed from
the alignment and the sequence. The number of soft clipped bases is added as
the 'ZA:i' (5' clipping) and 'ZB:i' (3' clipping) tags. The percentage of soft
clipped bases is given in the 'ZC:f' tag.
'''
import sys
import os
import pysam

def bam_iter(bam, quiet=True, show_ref_pos=False, ref=None, start=None, end=None, callback=None):
    '''
    >>> [x.qname for x in bam_iter(bam_open(os.path.join(os.path.dirname(__file__), 't', 'test.bam')), quiet=True)]
    ['A', 'B', 'E', 'C', 'D', 'F', 'Z']
    '''

    if os.path.exists('%s.bai' % bam.filename):
        # This is an indexed file, so it is ref sorted...
        # Meaning that we should show chrom:pos, instead of read names
        show_ref_pos = True

    eta = None

    if not ref:
        if not quiet and bam.filename:
            eta = ETA(os.stat(bam.filename).st_size)

        for read in bam:
            pos = bam.tell()
            bgz_offset = pos >> 16

            if not quiet and eta:
                if callback:
                    eta.print_status(bgz_offset, extra=callback(read))
                elif (show_ref_pos):
                    if read.tid > -1:
                        eta.print_status(bgz_offset, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))
                    else:
                        eta.print_status(bgz_offset, extra='unmapped %s' % (read.qname))
                else:
                    eta.print_status(bgz_offset, extra='%s' % read.qname)

            yield read

    else:
        working_chrom = None
        if ref in bam.references:
            working_chrom = ref
        elif ref[0:3] == 'chr':
            # compensate for Ensembl vs UCSC ref naming
            if ref[3:] in bam.references:
                working_chrom = ref[3:]

        if not working_chrom:
            raise ValueError('Missing reference: %s' % ref)

        tid = bam.gettid(working_chrom)

        if not start:
            start = 0
        if not end:
            end = bam.lengths[tid]

        if not quiet and bam.filename:
            eta = ETA(end - start)

        for read in bam.fetch(working_chrom, start, end):
            if not quiet and eta:
                if callback:
                    eta.print_status(read.pos - start, extra=callback(read))
                else:
                    eta.print_status(read.pos - start, extra='%s:%s %s' % (bam.getrname(read.tid), read.pos, read.qname))

            yield read

    if eta:
        eta.done()

def bam_removeclipping(infile, outfile):
    bam = pysam.Samfile(infile, "rb")
    out = pysam.Samfile(outfile, "wb", template=bam)
    total = 0
    count = 0
    unmapped = 0
    for read in bam_iter(bam):
        code = read_removeclipping(read)

        if code == 1:
            unmapped += 1
        elif code == 2:
            count += 1

        total += 1
        out.write(read)

    bam.close()
    out.close()
    sys.stderr.write('Wrote: %s reads\nAltered: %s\nUnmapped: %s\n' % (total, count, unmapped))


def read_removeclipping(read):
    if read.is_unmapped:
        return 1

    newcigar = []
    clip_5 = 0
    clip_3 = 0

    changed = False
    inseq = False
    for op, length in read.cigar:
        if op == 5:  # H
            changed = True
        elif op == 4:  # S
            changed = True
            if not inseq:
                clip_5 = length
            else:
                clip_3 = length
        else:
            inseq = True
            newcigar.append((op, length))

    if not changed:
        return 0

    read.cigar = newcigar
    orig_length = len(read.seq)

    s = read.seq
    q = read.qual

    if clip_3:
        read.seq = s[clip_5:-clip_3]
        if q:
            read.qual = q[clip_5:-clip_3]
    else:
        read.seq = s[clip_5:]
        if q:
            read.qual = q[clip_5:]

    newtags = []
    if clip_5:
        newtags.append(('ZA', clip_5))
    if clip_3:
        newtags.append(('ZB', clip_3))

    newtags.append(('ZC', float(clip_5 + clip_3) / orig_length))

    read.tags = read.tags + newtags

    return 2


def usage():
    print(__doc__)
    print("Usage: bamutils removeclipping inbamfile outbamfile")
    sys.exit(-1)

if __name__ == "__main__":
    infile = None
    outfile = None
    force = False

    for arg in sys.argv[1:]:
        if arg == "-h":
            usage()
#         elif arg == "-f":
#             force = True
        elif not infile:
            if os.path.exists(os.path.expanduser(arg)):
                infile = os.path.expanduser(arg)
            else:
                sys.stderr.write("File: %s not found!" % arg)
                usage()
        elif not outfile:
            if force or not os.path.exists(os.path.expanduser(arg)):
                outfile = arg
            else:
                sys.stderr.write(
                    "File: %s exists! Not overwriting without -f force." % arg
                    )
                usage()
        else:
            usage()

    if not infile or not outfile:
        usage()

    bam_removeclipping(infile, outfile)
