from ._utils import *
from rich.progress import track
from .consensus_variants import *
class Circle(object):
    def __init__(self, chr='', start=0, end=0, strand='', coverage=0, circ='', rawfrag=1, unit=(), readname='',
                 cigar='', unitSite=()):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand
        self.coverage = coverage
        self.circle = circ
        self.finalCirc = circ
        self.rawfrag = rawfrag
        self.unit = [unit]
        self.readname = [readname]
        self.cigar = [cigar]
        self.unitSite = [unitSite]
        self.consensus = consensusSequence()

    def isIdentical(self, circ):
        '''
        Checks if a candidate eccDNA is identical to the current eccDNA and merges relevant information if they are identical.
        '''
        if self.circle == circ.circle:
            if circ.readname[0] not in self.readname:
                self.coverage += circ.coverage
            self.readname.extend(circ.readname)
            self.unit.extend(circ.unit)
            self.cigar.extend(circ.cigar)
            self.unitSite.extend(circ.unitSite)
            return True
        else:
            return False

    def printf(self):
        return f'{self.chr}:{self.start}-{self.end}({self.strand})', f'{self.end - self.start + 1}', \
            f'{self.chr}_{self.start}_{self.end}_{self.strand}'

    def valueUnit(self, unit, cigar, readname, unitSite):
        self.unit = [unit]
        self.cigar = [cigar]
        self.readname = [readname]
        self.unitSite = [unitSite]

    def getConsensus(self):
        refSite = (self.chr, self.start, self.end, self.strand)
        self.consensus = self.consensus.getConsensus(ref=refSite, readsID=self.readname, \
                                                     unit=self.unit, cigar=self.cigar, readsSite=self.unitSite)
        seq, var = self.consensus.getSeqVar(refSite, self.finalCirc)
        return seq, var, len(seq)

    def initialize(self):
        self.unit = []
        self.cigar = []
        self.readname = []
        self.unitSite = []
        return self

class MultiCircle(object):
    def __init__(self):
        self.fragments = np.array([])
        self.circle = ''
        self.raw = []
        self.rawCoverage = []
        self.coverage = 0
        self.length = 0
        self.readname = []
        self.rawIdx = []

    def splitCircle(self, circle, raw, cov, nfrag=1, unit=(), readname='', indices='', cigar='', unitSite=()):
        self.circle = circle
        self.raw.append(raw)
        self.coverage = cov
        self.rawCoverage.append(cov)
        self.readname.append([readname])
        idx = list(map(int, indices.split('|')))
        for i, frag in enumerate(circle.split('|')):
            self.rawIdx.append([])
            tmp = parse_genomic_location(frag)
            self.length += tmp[2] - tmp[1] + 1
            self.fragments = np.append(self.fragments, Circle(chr=tmp[0], start=tmp[1], end=tmp[2], strand=tmp[3],
                                                              coverage=cov, circ=frag, rawfrag=nfrag + 1))
        idxF = idx.index(nfrag)
        self.fragments[idxF].valueUnit(unit, cigar, readname, unitSite)
        self.rawIdx[idxF] = [1]
        return self

    def isIdentical(self, circ):
        '''
        Checks if a multi-fragments candidate eccDNA is identical to the current eccDNA and merges relevant information if they are identical.
        '''
        if self.circle == circ.circle:
            if circ.raw[0] not in self.raw:
                self.raw.extend(circ.raw)
                self.rawCoverage.append(circ.coverage)
                self.readname.extend(circ.readname)
                self.coverage += circ.coverage
                rawIdx = len(self.readname)
            else:
                idx_f = self.raw.index(circ.raw[0])
                rawIdx = idx_f + 1
                if circ.readname[0][0] not in self.readname[idx_f]:
                    self.rawCoverage[idx_f] += circ.coverage
                    self.coverage += circ.coverage
                    self.readname[idx_f].extend(circ.readname[0])
            self.rawIdx[circ.rawIdx.index([1])].append(rawIdx)
            for i, single in enumerate(circ.fragments):
                self.fragments[i].unit.extend(single.unit)
                self.fragments[i].unit = list(filter(lambda x: x != (), self.fragments[i].unit))
                self.fragments[i].cigar.extend(single.cigar)
                self.fragments[i].cigar = list(filter(lambda x: x != '', self.fragments[i].cigar))
                self.fragments[i].readname.extend(single.readname)
                self.fragments[i].readname = list(filter(lambda x: x != '', self.fragments[i].readname))
                self.fragments[i].unitSite.extend(single.unitSite)
                self.fragments[i].unitSite = list(filter(lambda x: x != (), self.fragments[i].unitSite))
            return True
        else:
            return False

    def printf(self):
        line = []
        length = 0
        seqName = []
        for i in range(len(self.fragments)):
            circle, l, seqID = self.fragments[i].printf()
            line.append(circle)
            length += int(l)
            seqName.append(seqID)
        return f"{'|'.join(line)}", f'{length}', f"{'/'.join(seqName)}"

    def printfRaw(self):
        for circle in self.raw:
            yield f'{circle},{self.length}'

    def getConsensus(self):
        sequence = []
        variants = []
        seqLength = 0
        for circ in self.fragments:
            seq, var, length = circ.getConsensus()
            seqLength += length
            sequence.append(seq)
            variants.extend(var)
        sequence = ''.join(sequence)
        return sequence, variants, seqLength

class Circgroup(object):
    def __init__(self):
        self.group = []
        self.circ = Circle()

    def isEmpty(self):
        return len(self.group) == 0

    def compare_position(self, i, newCirc):
        '''
        Checks if the specified eccDNA in the current group matches the new eccDNA within the allowed fluctuate.
        '''
        return (self.group[i].chr == newCirc.chr) and \
            ((abs(self.group[i].start - newCirc.start) <= args.fluctuate) or
             (abs(self.group[i].end - newCirc.end) <= args.fluctuate))

    def multiCompare(self, i, newCirc):
        '''
        Checks if the specified multi-fragment eccDNA in the current group matches the new multi-fragment eccDNA within the allowed fluctuate.
        '''
        if self.group[0].fragments[0].strand == newCirc.fragments[0].strand:
            for j in range(len(newCirc.fragments)):
                if (self.group[i].fragments[j].chr == newCirc.fragments[j].chr) and \
                        ((abs(self.group[i].fragments[j].start - newCirc.fragments[j].start) <= args.fluctuate) or
                         (abs(self.group[i].fragments[j].end - newCirc.fragments[j].end) <= args.fluctuate)):
                    continue
                else:
                    return False
            return True
        else:
            for j in range(len(newCirc.fragments)):
                if (self.group[i].fragments[j].chr == newCirc.fragments[-j].chr) and \
                        ((abs(self.group[i].fragments[j].start - newCirc.fragments[-j].start) <= args.fluctuate) or
                         (abs(self.group[i].fragments[j].end - newCirc.fragments[-j].end) <= args.fluctuate)):
                    continue
                else:
                    return False
            return True

    def inclusive(self, i, newCirc):
        return (self.group[i].chr == newCirc.chr) and \
            ((abs(self.group[i].start - newCirc.start) <= args.fluctuate) and
             (abs(self.group[i].end - newCirc.end) <= args.fluctuate))

    def judgeRegion(self, circle, single=True):
        if self.isEmpty():
            return True
        if single:
            res = [self.compare_position(i, circle) for i in range(len(self.group))]
            if True in res:
                return True
            else:
                return False
        else:
            if len(circle.fragments) != len(self.group[0].fragments):
                return False
            res = [self.multiCompare(i, circle) for i in range(len(self.group))]
            if True in res:
                return True
            else:
                return False

    def multiInclusive(self, i, newCirc):
        for j in range(len(newCirc.fragments)):
            if (self.group[i].fragments[j].chr == newCirc.fragments[j].chr) and \
                    (self.group[i].fragments[j].strand == newCirc.fragments[j].strand) and \
                    ((abs(self.group[i].fragments[j].start - newCirc.fragments[j].start) <= args.fluctuate) and \
                     (abs(self.group[i].fragments[j].end - newCirc.fragments[j].end) <= args.fluctuate)):
                continue
            else:
                return False
        return True

    def judgeGroup(self, circle, single=True):
        if self.isEmpty():
            return True
        if single:
            res = [self.inclusive(i, circle) for i in range(len(self.group))]
        else:
            res = [self.multiInclusive(i, circle) for i in range(len(self.group))]
        if True in res:
            return True
        else:
            return False

    def find_longest_list(self, groups):
        longest_group = []
        longest = 0
        for group in groups:
            if len(group) > longest:
                longest = len(group)
                longest_group = group
        return longest_group

    def remove_extremum(self, data):
        data.sort()
        groups = []
        i = 0
        n = len(data)
        while i < n:
            current = data[i]
            in_group = False
            for m in range(len(groups)):
                group = groups[m]
                if min(group) - 20 <= current and max(group) + 20 >= current:
                    in_group = True
                    group.append(current)
                    groups[m] = group
                    break
            if not in_group:
                new_group = [current]
                low = current
                high = current
                j = i + 1
                while j < n and low - 20 <= data[j] and high + 20 >= data[j]:
                    high = max(high, data[j])
                    low = min(low, data[j])
                    new_group.append(data[j])
                    j += 1
                i = j
                groups.append(new_group)
        return self.find_longest_list(groups)

    def single_fragment_ecc(self):
        self.circ.coverage += sum([group.coverage for group in self.group])
        if self.circ.coverage >= args.cov and len(self.group) >= args.nc:
            self.circ.unit = [unit for group in self.group for unit in group.unit]
            self.circ.readname = [readname for group in self.group for readname in group.readname]
            self.circ.cigar = [cigar for group in self.group for cigar in group.cigar]
            self.circ.unitSite = [unitSite for group in self.group for unitSite in group.unitSite]
            unit = [u for group in self.group for u in group.unit]
            self.circ.chr = unit[0][0]
            self.circ.strand = unit[0][3]
            self.circ.start = min(multimode(self.remove_extremum(np.array(unit)[:, 1].astype(int))))
            self.circ.end = max(multimode(self.remove_extremum(np.array(unit)[:, 2].astype(int))))
            self.circ.finalCirc = f'{self.circ.chr}:{self.circ.start}-{self.circ.end}({self.circ.strand})'
            return True
        else:
            return False

    def multi_fragments_ecc(self):
        self.circ = MultiCircle()
        self.circ.coverage += sum([group.coverage for group in self.group])
        self.circ.readname = [r for group in self.group for readname in group.readname for r in readname]
        if self.circ.coverage >= args.cov and (len(self.group) >= args.nc or len(
                np.unique([r for group in self.group for r in group.raw])) >= args.nc):
            f = []
            for i, _ in enumerate(self.group[0].fragments):
                unit = [u for group in self.group for u in group.fragments[i].unit]
                chr_ = unit[0][0]
                strand_ = unit[0][3]
                start_ = min(multimode(self.remove_extremum(np.array(unit)[:, 1].astype(int))))
                end_ = max(multimode(self.remove_extremum(np.array(unit)[:, 2].astype(int))))
                f.append(f'{chr_}:{start_}-{end_}({strand_})')
                self.circ.fragments = np.append(self.circ.fragments,
                                                Circle(chr=chr_, start=start_, end=end_, strand=strand_).initialize())
            f = '|'.join(f)
            for group in self.group:
                for i in range(len(group.fragments)):
                    self.circ.fragments[i].finalCirc = f
                    self.circ.fragments[i].cigar.extend(group.fragments[i].cigar)
                    self.circ.fragments[i].unit.extend(group.fragments[i].unit)
                    self.circ.fragments[i].unitSite.extend(group.fragments[i].unitSite)
                    self.circ.fragments[i].readname.extend(group.fragments[i].readname)
                    self.circ.fragments[i].chr = group.fragments[i].chr
            return True
        else:
            return False

class Region(object):
    def __init__(self, frag, raw, depth, chr_, start, end, strand, reads, cigars, unitStart, unitEnd):
        self.frag = frag
        self.rawFrag = raw
        self.depth = depth
        self.chr_ = chr_
        self.start = start
        self.end = end
        self.strand = strand
        self.reads = reads
        self.cigar = cigars
        self.unitStart = unitStart
        self.unitEnd = unitEnd

    def single_fragment(self):
        lastCirc = Circle()
        plusCirc = Circgroup()
        minusCirc = Circgroup()
        for circ, cov, chr_, start_, end_, strand_, readname, cigar_, us, ue in \
                track(zip(self.frag, self.depth, self.chr_, self.start, self.end, self.strand, \
                          self.reads, self.cigar, self.unitStart, self.unitEnd),
                      description='Get accurate location from single fragment candidate eccDNA...'):
            tmp = parse_genomic_location(circ)
            single = Circle(chr=tmp[0], start=tmp[1], end=tmp[2], strand=tmp[3], coverage=cov, circ=circ,
                            unit=(chr_, start_, end_, strand_), readname=readname, cigar=cigar_, unitSite=(us, ue))
            if lastCirc.isIdentical(single):
                continue
            lastCirc = single
            if single.strand == '+':
                if plusCirc.judgeRegion(single):
                    plusCirc.group.append(single)
                else:
                    yield plusCirc
                    plusCirc = Circgroup()
                    plusCirc.group.append(single)
            else:
                if minusCirc.judgeRegion(single):
                    minusCirc.group.append(single)
                else:
                    yield minusCirc
                    minusCirc = Circgroup()
                    minusCirc.group.append(single)
        if not plusCirc.isEmpty():
            yield plusCirc
        if not minusCirc.isEmpty():
            yield minusCirc

    def multi_fragments(self, Nfragment, indices):
        lastCirc = MultiCircle()
        circGroup = Circgroup()
        for circ, raw, cov, chr_, start_, end_, strand_, readname, circIdx, fragIdx, cigar_, us, ue in \
                track(zip(self.frag, self.rawFrag, self.depth, self.chr_, self.start, self.end, self.strand, \
                          self.reads, indices, Nfragment, self.cigar, self.unitStart, self.unitEnd), \
                      description='Get accurate locations from multiple fragments candidate eccDNA...'):
            multi = MultiCircle().splitCircle(circ, raw, cov, nfrag=fragIdx, unit=(chr_, start_, end_, strand_), \
                                              readname=readname, indices=circIdx, cigar=cigar_, unitSite=(us, ue))
            if lastCirc.isIdentical(multi):
                continue
            lastCirc = multi
            if circGroup.judgeRegion(multi, single=False):
                circGroup.group.append(multi)
            else:
                yield circGroup
                circGroup = Circgroup()
                circGroup.group.append(multi)
        if not circGroup.isEmpty():
            yield circGroup

class SplitRegion(object):
    def __init__(self, region):
        self.region = region
        self.group = [Circgroup()]
        self.printfInfo = []
        self.sequence = []
        self.variants = []

    def merge(self, idx):
        for i in idx[1:]:
            self.group[idx[0]].group.extend(self.group[i].group)
        self.group[idx[0]].group = list(set(self.group[idx[0]].group))
        self.group = [element for i, element in enumerate(self.group) if i not in idx[1:]]

    def divide(self, single=True):
        for i in range(len(self.region.group)):
            idx = []
            for j in range(len(self.group)):
                circ = self.group[j]
                if circ.judgeGroup(self.region.group[i], single=single):
                    idx.append(j)
                    circ.group.append(self.region.group[i])
            if len(idx) == 0:
                tmp = Circgroup()
                tmp.group.append(self.region.group[i])
                self.group.append(tmp)
            elif len(idx) > 1:
                self.merge(idx)
        return self

    def judge(self):
        if (np.array([x.start for x in self.region.group]).ptp() <= args.fluctuate) and \
                (np.array([x.end for x in self.region.group]).ptp() <= args.fluctuate):
            self.group = [self.region]
            return self
        else:
            return self.divide()

    def single_printf(self):
        for group in self.group:
            if group.single_fragment_ecc():
                circ_info = self._process_group_circ(group)
                yield from self._process_consolidate_single_candidate(group, circ_info)
            else:
                yield from self._process_non_consolidate_single(group.group)

    def multi_printf(self):
        for group in self.group:
            if group.multi_fragments_ecc():
                circ_info = self._process_group_circ(group)
                yield from self._process_consolidate_multi_candidate(group, circ_info)
            else:
                yield from self._process_non_consolidate_multi(group.group)

    def _process_group_circ(self, group):
        sequence, variants, seqLength = group.circ.getConsensus()
        self.variants.extend(variants)
        circSite, length, seqName = group.circ.printf()
        self.sequence.append(f'>{seqName}\n{sequence}')
        return (circSite, length, seqLength)

    def _process_single_circ(self, circ):
        sequence, variants, seqLength = circ.getConsensus()
        self.variants.extend(variants)
        circSite, length, seqName = circ.printf()
        self.sequence.append(f'>{seqName}\n{sequence}')
        return (circSite, length, seqLength)

    def _append_printf_info(self, circSite, cov, reads, length, seqLength):
        unique_reads = len(np.unique(reads))
        count = circSite.count("|") + 1
        self.printfInfo.append(
            f"{circSite},{cov},{count},{unique_reads},{length},{seqLength}"
        )

    def _process_consolidate_single_candidate(self, group, circ_info):
        circSite, length, seqLength = circ_info
        reads = []
        cov = 0
        for circ in group.group:
            reads.extend(circ.readname)
            cov += circ.coverage
            rawSite, rawLen, _ = circ.printf()
            yield f"{circSite},{length},True,{rawSite},{rawLen},{circ.coverage}\n"
        self._append_printf_info(circSite, cov, reads, length, seqLength)

    def _process_non_consolidate_single(self, circs):
        for circ in circs:
            circSite, length, seqLength = self._process_single_circ(circ)
            self._append_printf_info(circSite, circ.coverage, circ.readname, length, seqLength)
            tmp = f"{circSite},{length}"
            yield f"{tmp},False,{tmp},{circ.coverage}\n"

    def _process_consolidate_multi_candidate(self, group, circ_info):
        circSite, length, seqLength = circ_info
        reads = []
        cov = 0
        for circ in group.group:
            reads.extend(r for sublist in circ.readname for r in sublist)
            cov += circ.coverage
            for n, line in enumerate(circ.printfRaw()):
                yield f"{circSite},{length},True,{line},{circ.rawCoverage[n]}\n"
        self._append_printf_info(circSite, cov, reads, length, seqLength)

    def _process_non_consolidate_multi(self, circs):
        for circ in circs:
            for i, raw in enumerate(circ.raw):
                tmpCirc = self._construct_tmp_circ(circ, i, raw)
                circ_info = self._process_single_circ(tmpCirc)
                circSite, length, seqLength = circ_info
                self._append_printf_info(circSite, tmpCirc.coverage, tmpCirc.readname, length, seqLength)
                line = next(tmpCirc.printfRaw())
                yield f"{circSite},{length},False,{line},{tmpCirc.rawCoverage[0]}\n"

    def _construct_tmp_circ(self, circ, i, raw):
        tmpCirc = MultiCircle()
        tmpCirc.circle = raw
        tmpCirc.raw.append(raw)
        tmpCirc.rawCoverage.append(circ.rawCoverage[i])
        tmpCirc.readname = circ.readname[i]
        tmpCirc.coverage = circ.rawCoverage[i]
        tmpCirc.length = circ.length
        tmp = raw.split('|')
        tmpCirc.fragments = np.array([Circle() for _ in range(len(tmp))])
        minIdx = np.argmin(tmp)
        tmpIdx = [*range(len(tmp))[minIdx:], *range(len(tmp))[:minIdx]]
        for j, idx in enumerate(tmpIdx):
            fragment = circ.fragments[idx]
            tmpCirc.fragments[j].finalCirc = raw
            tmpCirc.fragments[j].chr = fragment.chr
            tmpCirc.fragments[j].start = fragment.start
            tmpCirc.fragments[j].end = fragment.end
            tmpCirc.fragments[j].strand = fragment.strand

            idxUnit = np.array(circ.rawIdx[idx]) == i + 1
            tmpCirc.fragments[j].cigar = np.array(fragment.cigar)[idxUnit].tolist()
            tmpCirc.fragments[j].unit = [tup for tup, bool_ in zip(fragment.unit, idxUnit) if bool_]
            tmpCirc.fragments[j].readname = np.array(fragment.readname)[idxUnit].tolist()
            tmpCirc.fragments[j].unitSite = np.array(fragment.unitSite)[idxUnit].tolist()
        return tmpCirc

    def writeInfo(self):
        return '\n'.join(self.printfInfo), '\n'.join(self.sequence)

    def getVariants(self):
        for var in self.variants:
            if var.strand == '-':
                var.ref = str(Bio.Seq.Seq(var.ref).reverse_complement())
                var.alt = str(Bio.Seq.Seq(var.alt).reverse_complement())
            var = var.varType()
            yield ','.join([var.chr, str(var.pos), var.ref, var.alt, \
                             str(var.count), str(var.depth), var.type, var.site]) + '\n'
