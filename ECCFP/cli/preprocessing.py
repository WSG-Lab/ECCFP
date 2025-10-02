import sys
from ._utils import *
from .config import args
from rich.progress import track

class Interval(object):
    def __init__(self, chr, start, end, strand):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = strand

class Fragment(object):
    def __init__(self, start, end, chr, chrStart, chrEnd, strand, mapQual, cigar, readname='', readUnit=False):
        self.start = start
        self.end = end
        self.interval = Interval(chr, chrStart, chrEnd, strand)
        self.mapQual = mapQual
        self.cigar = cigar
        self.name = readname
        self.readUnit = readUnit
        self.isGhost = False
        self.number = 0

class FragmentGroup(object):
    def __init__(self, read='', seqname=''):
        self.group = []
        self.intervals = []
        self.nfrag = 0
        self.passes = 0
        self.pointer = 0
        self.addFirst = False
        self.addLast = False
        self.addpass = False
        self.maxOffset = args.maxOffset
        self.readname = read
        self.seqname = seqname

    def _sameLocation(self, index, interval, withStart=True, withEnd=True, Inclusive=False, terminal=False):
        """
            Determines if the given fragment is at the same location as the eccDNA fragment.
        """
        if terminal:
            return False
        if len(self.group[index]) == 0:
            return True
        if interval.chr != self.group[index][0].interval.chr:
            return False
        if interval.strand != self.group[index][0].interval.strand:  ##### question1
            return False
        if withStart:
            start = multimode([frag.interval.start for frag in self.group[index]])
            res = [abs(i - interval.start) > self.maxOffset for i in start]
            if False not in res:
                return False
        if withEnd:
            end = multimode([frag.interval.end for frag in self.group[index]])
            res = [abs(i - interval.end) > self.maxOffset for i in end]
            if False not in res:
                return False
        if Inclusive:
            start = multimode([frag.interval.start for frag in self.group[index]])
            end = multimode([frag.interval.end for frag in self.group[index]])
            if interval.start < min(start) - args.maxOffset or interval.end > max(end) + args.maxOffset:
                return False
        return True

    def _rewind(self, interval):
        return self._sameLocation(0, interval)

    def _noOverlap(self, fragment):
        interval = fragment.interval
        n = 0
        for frag in self.group:
            if frag[0].interval.chr == interval.chr:
                if abs(frag[-1].end - fragment.start) > 1:
                    n += 1
                if interval.start <= frag[0].interval.end and interval.end >= frag[0].interval.start and \
                        interval.strand == frag[0].interval.strand:
                    return False
            else:
                if abs(frag[-1].end - fragment.start) > 1:
                    n += 1
        if n == len(self.group):
            return False
        return True

    def add(self, fragment):
        if len(self.group) == 0:
            self.group.append([fragment])
            return True
        if self.nfrag == 0:
            if self._rewind(fragment.interval):
                self.group[0].append(fragment)
                self.nfrag = len(self.group)
                self.passes = 1
                if self.pointer == 0:
                    self.passes += 1
                self.pointer = 0
                return True
            else:
                if self._noOverlap(fragment):
                    self.group.append([fragment])
                    self.pointer += 1
                    return True
        else:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
            if self._sameLocation(pointer, fragment.interval):
                self.group[pointer].append(fragment)
                self.pointer = pointer
                if self.pointer == self.nfrag - 1:
                    self.passes += 1
                return True
        return False

    def addTerminals(self, firstFragment, lastFragment, first, last):
        start = min(multimode([frag.interval.start for frag in self.group[-1]]))
        end = max(multimode([frag.interval.end for frag in self.group[-1]]))
        if len(self.group) == 0:
            return False, False
        self.addLast = False
        self.addFirst = False
        # Add last fragment
        if self.nfrag > 0:
            pointer = self.pointer + 1
            if pointer == self.nfrag:
                pointer = 0
        else:
            pointer = 0
        if lastFragment.interval.strand == '+':
            withStart = True
            withEnd = False
        else:
            withStart = False
            withEnd = True
        if self._sameLocation(pointer, lastFragment.interval, withStart, withEnd, Inclusive=True, terminal=last):
            self.addLast = True
            self.pointer = pointer
        # Add first fragment
        if firstFragment.interval.strand == '+':
            withStart = False
            withEnd = True
        else:
            withStart = True
            withEnd = False
        if self._sameLocation(-1, firstFragment.interval, withStart, withEnd, Inclusive=True, terminal=first):
            self.addFirst = True
        if self.addLast and not last:
            self.group[pointer].append(lastFragment)
            last = True
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        else:
            self.addLast = False
        if self.addFirst and not first:
            self.group[-1].append(firstFragment)
            first = True
            if self.nfrag == 0:
                self.nfrag = len(self.group)
                self.passes = 1
        else:
            self.addFirst = False

        # Determine whether the first fragment and the last fragment can form a complete full pass
        if self.addFirst and self.nfrag == 1:
            if abs(firstFragment.interval.start - start) <= 20 and abs(firstFragment.interval.end - end) <= 20:
                self.passes += 1
                self.addpass = True
        if self.addLast and self.nfrag == 1:
            if pointer == self.nfrag - 1:
                if abs(lastFragment.interval.start - start) <= 20 and abs(lastFragment.interval.end - end) <= 20:
                    self.passes += 1
                    self.addpass = True
        if self.addLast and self.addFirst and not self.addpass:
            if pointer == self.nfrag - 1 and firstFragment.interval.strand == lastFragment.interval.strand:
                if firstFragment.interval.strand == '+':
                    if firstFragment.interval.start <= lastFragment.interval.end + 1:
                        self.addpass = True
                        self.passes += 1
                else:
                    if lastFragment.interval.start <= firstFragment.interval.end + 1:
                        self.addpass = True
                        self.passes += 1
        return first, last

    def extract_position(self):
        if self.nfrag == 0:
            return None
        for frags in self.group:
            chr = frags[0].interval.chr
            start = min(multimode([frag.interval.start for frag in frags]))
            end = max(multimode([frag.interval.end for frag in frags]))
            strand = multimode([frag.interval.strand for frag in frags])[0]
            self.intervals.append(Interval(chr, start, end, strand))
        return self

    def merge_candidate(self, *fraggroup):
        '''
        Merges candidate eccDNA fragments from the same position but non-consecutive regions
        to form the final candidate eccDNA.
        '''
        self.nfrag = fraggroup[0].nfrag
        l = lambda x: x.nfrag == self.nfrag
        nfragjud = [l(i) for i in fraggroup]
        if False in nfragjud:
            sys.stdout.write('Error!')
            sys.exit()
        for i in range(self.nfrag):
            self.group.append([g for frag in fraggroup for g in frag.group[i]])
            self.group[i].sort(key=lambda x: x.start)
        self.passes = sum([i.passes for i in fraggroup])
        for frag in fraggroup:
            if frag.passes:
                self.addpass = True
            if frag.addFirst and not self.addFirst:
                self.addFirst = True
            if frag.addLast and not self.addLast:
                self.addLast = True
            if frag.addLast and frag.addFirst:
                break
        if self.addFirst and self.addLast and not self.addpass and len(self.group) == 1:
            firstFragment = self.group[0][0]
            lastFragment = self.group[-1][-1]
            if firstFragment.interval.strand == lastFragment.interval.strand:
                if firstFragment.interval.strand == '+':
                    if firstFragment.interval.start <= lastFragment.interval.end + 1:
                        self.passes += 1
                else:
                    if lastFragment.interval.start <= firstFragment.interval.end + 1:
                        self.passes += 1
        return self.extract_position()

    def adjust_fragments(self, index):
        self.group = list(np.array(self.group, dtype=object)[index])
        self.intervals = list(np.array(self.intervals, dtype=object)[index])
        return self

    def output_unit_info(self, read):
        '''
        Outputs the information for unit.txt, including fragment details.
        '''
        frags = []
        for interval in self.intervals:
            frags.append(interval.chr + ':' + str(interval.start + 1) + '-' + str(
                interval.end + 1) + '(' + interval.strand + ')')
        frags = '|'.join(frags)
        n = 0
        unit = ''
        for group in self.group:
            n += 1
            for fragment in group:
                unit = f'{unit}{read.name}\t{self.passes}\t{n}\t{fragment.start + 1}\t{fragment.end + 1}\t' \
                       f'{fragment.interval.chr}\t{fragment.interval.start + 1}\t{fragment.interval.end + 1}\t' \
                        f'{fragment.interval.strand}\t{frags}\t{fragment.isGhost}\t{fragment.number}\t{fragment.cigar}\n'
        return unit

class Read(object):
    def __init__(self, name=None, length=0):
        self.name = name
        self.length = int(length)
        self.fragments = []
        self.minMapQual = args.minMapQual
        self.circle = []
        self.delFirst = False
        self.delLast = False
        self.eccSite = ''

    def isEmpty(self):
        return len(self.fragments) == 0

    def add(self, fragment):
        self.fragments.append(fragment)

    def ordered(self):
        self.fragments.sort(key=lambda x: x.start)
        return self

    def threadFix(self):
        '''
        Trims the overlapping parts of fragments.
        '''
        self.ordered()
        n = 0
        self.fragments[-1].number = len(self.fragments)
        for i in range(len(self.fragments) - 1):
            n += 1
            self.fragments[i].number = n
            dislen = self.fragments[i + 1].start - self.fragments[i].end - 1
            if dislen > 0:
                continue
            elif dislen < 0:
                dislen = -dislen
                if self.fragments[i].interval.strand == '+':
                    self.fragments[i].end -= dislen
                    cigar_re = re.compile(r"(\d+)([MDI])")
                    ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i].cigar)]
                    p = len(ces) - 1
                    while dislen > 0 and p >= 0:
                        if ces[p][1] == 'M' or ces[p][1] == 'I':
                            s = min(ces[p][0], dislen)
                            dislen -= s
                            ces[p][0] -= s
                            if ces[p][1] == 'M':
                                self.fragments[i].interval.end -= s
                            if ces[p][0] == 0:
                                p -= 1
                        if p >= 0 and ces[p][1] == 'D':
                            self.fragments[i].interval.end -= ces[p][0]
                            ces[p][0] = 0
                            p -= 1
                    if p < 0:
                        return False
                    if ces[p][1] == 'I':
                        ces[p][1] = 'M'
                        self.fragments[i].interval.end += ces[p][0]
                    self.fragments[i].cigar = ''.join([str(step) + op for step, op in ces[0:p + 1]])
                else:
                    if self.fragments[i + 1].interval.strand == '-':
                        self.fragments[i + 1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i + 1].cigar)]
                        p = len(ces) - 1
                        while dislen > 0 and p >= 0:
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i + 1].interval.end -= s
                                if ces[p][0] == 0:
                                    p -= 1
                            if p >= 0 and ces[p][1] == 'D':
                                self.fragments[i + 1].interval.end -= ces[p][0]
                                ces[p][0] = 0
                                p -= 1
                        if p < 0:
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i + 1].interval.end += ces[p][0]
                        self.fragments[i + 1].cigar = ''.join([str(step) + op for step, op in ces[0:p + 1]])
                    else:
                        self.fragments[i + 1].start += dislen
                        cigar_re = re.compile(r"(\d+)([MDI])")
                        ces = [[int(m.group(1)), m.group(2)] for m in cigar_re.finditer(self.fragments[i + 1].cigar)]
                        p = 0
                        while dislen > 0 and p < len(ces):
                            if ces[p][1] == 'M' or ces[p][1] == 'I':
                                s = min(ces[p][0], dislen)
                                dislen -= s
                                ces[p][0] -= s
                                if ces[p][1] == 'M':
                                    self.fragments[i + 1].interval.start += s
                                if ces[p][0] == 0:
                                    p += 1
                            if p < len(ces) and ces[p][1] == 'D':
                                self.fragments[i + 1].interval.start += ces[p][0]
                                ces[p][0] = 0
                                p += 1
                        if p >= len(ces):
                            return False
                        if ces[p][1] == 'I':
                            ces[p][1] = 'M'
                            self.fragments[i + 1].interval.start -= ces[p][0]
                        self.fragments[i + 1].cigar = ''.join([str(step) + op for step, op in ces[p:]])
        return True

    def remove_contained_fragments(self):
        '''
        Removes shorter fragments that are completely contained within longer fragments.

        This function iterates through the list of fragments and
        checks if one fragment is completely contained within another.
        The shorter fragment is marked for deletion.
        '''
        deletion = []
        for i in range(len(self.fragments) - 1):
            if self.fragments[i].start >= self.fragments[i + 1].start and \
                    self.fragments[i].end <= self.fragments[i + 1].end:
                deletion.append(self.fragments[i])
            elif self.fragments[i].start <= self.fragments[i + 1].start and \
                    self.fragments[i].end >= self.fragments[i + 1].end:
                deletion.append(self.fragments[i + 1])
            else:
                continue
        if self.fragments[0] in deletion:
            self.delFirst = True
        if self.fragments[-1] in deletion:
            self.delLast = True
        if deletion:
            self.fragments = [frag for frag in self.fragments if frag not in deletion]
        return self.threadFix()

    def compare_candidate(self, newCirc, lastCirc, lenF):
        for n in range(lenF):
            newLoc = parse_genomic_location(newCirc[n])
            lastLoc = parse_genomic_location(lastCirc[n])
            if newLoc[0] == lastLoc[0] and newLoc[3] == lastLoc[3] and \
                    abs(newLoc[1] - lastLoc[1]) <= args.fluctuate and \
                    abs(newLoc[2] - lastLoc[2]) <= args.fluctuate:
                continue
            else:
                return False
        return True

    def get_candidate_group(self, sortCirc, rawCirc, indices):
        candidateGroup = []
        for i in range(len(sortCirc)):
            circ = sortCirc[i]
            fragments = tuple(circ.split('|'))
            lenF = len(fragments)
            if len(candidateGroup) == 0:
                candidateGroup.append((circ, lenF, rawCirc[i], indices[i]))
                continue
            if lenF != candidateGroup[0][1]:
                yield candidateGroup
                candidateGroup = [(circ, lenF, rawCirc[i], indices[i])]
                continue
            judgeSame = False
            for groupF in candidateGroup:
                groupC = groupF[0].split("|")
                if self.compare_candidate(fragments, groupC, lenF):
                    judgeSame = True
                    break
            if judgeSame:
                candidateGroup.append((circ, lenF, rawCirc[i], indices[i]))
            else:
                yield candidateGroup
                candidateGroup = [(circ, lenF, rawCirc[i], indices[i])]
        if len(candidateGroup) != 0:
            yield candidateGroup

    def check_offset(self, dic):
        '''
        Check and Group Candidate eccDNA by Offset
        '''
        sortCirc, rawCirc, indices = np.vectorize(fragment_priority_adjuster)(list(dic.keys()))
        lexSortIdx = np.lexsort((sortCirc, np.vectorize(len)(sortCirc)))
        sortCirc = sortCirc[lexSortIdx]
        rawCirc = rawCirc[lexSortIdx]
        indices = indices[lexSortIdx]
        group = []
        for circGroup in self.get_candidate_group(sortCirc, rawCirc, indices):
            tmp = []
            for circ, _, raw, idx in circGroup:
                tmp.extend([d.adjust_fragments(np.array(list(map(int, idx.split('|')))) - 1) for d in dic[raw]])
            same = False
            for g in group:
                if len(tmp[0].intervals) != len(g[0].intervals):
                    continue
                res = all([True if tmp[0].intervals[i].chr == g[0].intervals[i].chr and \
                                   tmp[0].intervals[i].strand == g[0].intervals[i].strand and \
                                   abs(tmp[0].intervals[i].start - g[0].intervals[i].start) <= args.fluctuate and \
                                   abs(tmp[0].intervals[i].end - g[0].intervals[i].end) <= args.fluctuate else False \
                           for i in range(len(g[0].intervals))])
                if res:
                    same = True
                    g.extend(tmp)
                    break
            if not same:
                group.append(tmp)
        return group

    def consolidate_candidate_groups(self, group):
        '''
        Consolidate candidate eccDNA groups.

        This function merges a list of FragmentGroup objects into fewer groups where each group contains FragmentGroup objects with similar locations.
        The consolidating process includes:
        1. Grouping FragmentGroup objects by their locations.
        2. Checking if the offset between FragmentGroup objects in each group is within the allowed range.
        3. Merging FragmentGroup objects with similar locations.
        '''
        dic = {}
        for frag in group:
            circle = []
            for interval in frag.intervals:
                circle.append(f'{interval.chr}:{interval.start}-{interval.end}({interval.strand[0]})')
            circle = '|'.join(circle)
            tmp = dic.get(circle, [])
            tmp.append(frag)
            dic[circle] = tmp
        if dic:
            group = self.check_offset(dic)
        else:
            group = []
        tmp = []
        for values in group:
            if len(values) == 1:
                tmp.append(values[0])
            else:
                tmp.append(FragmentGroup().merge_candidate(*values))
        return tmp

    def add_first_last(self, circle):
        tmp = []
        for fraggroup in circle:
            if self.delFirst:
                first = True
            else:
                for groupslice in fraggroup.group:
                    if self.fragments[1] in groupslice and \
                            self.fragments[1].start - self.fragments[0].end <= args.maxOffset:
                        first = False
                        break
                    else:
                        first = True
            if self.delLast:
                last = True
            else:
                for groupslice in fraggroup.group:
                    if self.fragments[-2] in groupslice and \
                            self.fragments[-1].start - self.fragments[-2].end <= args.maxOffset:
                        last = False
                        break
                    else:
                        last = True
            _, _ = fraggroup.addTerminals(self.fragments[0], self.fragments[-1], first, last)
            fraggroup = fraggroup.extract_position()
            if fraggroup:
                tmp.append(fraggroup)
        return self.consolidate_candidate_groups(tmp)

    def FragmentIsGhost(self, i):
        for j in range(i, len(self.fragments)):
            self.fragments[j].isGhost = True
        return True

    def bootstrap(self):
        if len(self.fragments) < 3:
            return None
        res = self.threadFix()
        if res == False:
            n = 0
            while not res:
                n += 1
                res = self.remove_contained_fragments()
                if (n > 2) or (len(self.fragments) < 3):
                    return None
        fraggroup = FragmentGroup()
        circle = []
        lastFrag = self.fragments[0]
        isGhost = False
        for i in range(1, len(self.fragments) - 1):
            if self.fragments[i].mapQual < self.minMapQual:
                if not isGhost:
                    isGhost = self.FragmentIsGhost(i)
                continue
            if self.fragments[i].start - lastFrag.end > args.maxOffset and len(fraggroup.group) != 0:
                res = False
            else:
                res = fraggroup.add(self.fragments[i])
            if res == False:
                if len(fraggroup.group) != 0:
                    if not isGhost:
                        isGhost = self.FragmentIsGhost(i)
                    circle.append(fraggroup)
                    fraggroup = FragmentGroup()
                    _ = fraggroup.add(self.fragments[i])
            lastFrag = self.fragments[i]
        if len(fraggroup.group) != 0:
            circle.append(fraggroup)
        self.circle = self.add_first_last(circle)
        return self.circle

class PAF(object):
    def __init__(self, filepath):
        self.filepath = filepath

    def reads(self):
        current = ''
        read = Read()
        with open(self.filepath, 'r') as paf:
            for line in track(paf, description='Indentify candidate eccDNA...'):
                fields = line.strip().split('\t')
                if fields[0] != current:
                    if not read.isEmpty():
                        yield read
                    current = fields[0]
                    read = Read(name=fields[0], length=int(fields[1]))
                assert fields[-1][0:2] == 'cg', 'The last column in PAF should be CIGAR string'
                frag = Fragment(start=int(fields[2]), end=int(fields[3]) - 1,
                                chr=fields[5], chrStart=int(fields[7]), chrEnd=int(fields[8]) - 1, strand=fields[4],
                                mapQual=int(fields[11]), cigar=fields[-1].split(':')[2])
                read.add(frag)
            if not read.isEmpty():
                yield read
