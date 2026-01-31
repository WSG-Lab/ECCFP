import os
from ECCFP.cli.accurate_eccDNA import *
from ECCFP.cli.preprocessing import PAF

def printf(unitPath, outLoc, outInfo, outSeq, outVar):
    a = pd.read_csv(unitPath)
    Nfullpass = a.Nfullpass.values
    idx = np.where(Nfullpass >= args.nf)
    Nfragment = a.fragments.values[idx]
    eccDNA = a.candidate_eccDNA.values[idx]
    Nfullpass = Nfullpass[idx]
    chr_ = a.chr.values[idx]
    start = a.start.values[idx]
    end = a.end.values[idx]
    strand = a.strand.values[idx]
    reads = a.reads.values[idx]
    eccDNA = eccDNA.astype(str)
    cigar = a.cigar.values[idx]
    read_start = a.read_start.values[idx]
    read_end = a.read_end.values[idx]
    with open(outLoc, 'w') as loc, open(outInfo, 'w') as out, open(outSeq, 'w') as fa, open(outVar, 'w') as var:
        loc.write(f'eccDNApos,eccDNA_len,consolidating,cand_eccDNA,cand_len,cand_Nfullpass\n')
        out.write(f'eccDNApos,Nfullpass,Nfragments,Nreads,refLength,seqLength\n')
        condation = np.char.find(eccDNA, '|') == -1
        idx = np.where(condation)
        frag = eccDNA[idx]
        if len(frag) > 0:
            cov = Nfullpass[idx]
            unit_chr_ = chr_[idx]
            unit_start = start[idx]
            unit_end = end[idx]
            unit_strand = strand[idx]
            unit_reads = reads[idx]
            unit_cigar = cigar[idx]
            unit_rs = read_start[idx]
            unit_re = read_end[idx]
            sortIndex = np.lexsort((frag, np.vectorize(len)(frag)))
            for region in Region(frag=frag[sortIndex], raw=frag[sortIndex], depth=cov[sortIndex], \
                                 chr_=unit_chr_[sortIndex], start=unit_start[sortIndex], \
                                 end=unit_end[sortIndex], strand=unit_strand[sortIndex], \
                                 reads=unit_reads[sortIndex], cigars=unit_cigar[sortIndex], \
                                 unitStart=unit_rs[sortIndex], unitEnd=unit_re[sortIndex]).single_fragment():
                result = SplitRegion(region).judge()
                for line in result.single_printf():
                    loc.write(line)
                info, seq = result.writeInfo()
                out.write(f'{info}\n')
                fa.write(f'{seq}\n')
                for v in result.getVariants():
                    var.write(v)
        idx = ~condation
        frag = eccDNA[idx]
        if len(frag) > 0:
            cov = Nfullpass[idx]
            unit_nfrag = Nfragment[idx]
            unit_chr_ = chr_[idx]
            unit_start = start[idx]
            unit_end = end[idx]
            unit_strand = strand[idx]
            unit_reads = reads[idx]
            unit_cigar = cigar[idx]
            unit_rs = read_start[idx]
            unit_re = read_end[idx]
            f = np.vectorize(fragment_priority_adjuster)
            frag = f(frag)
            sortIndex = np.lexsort((frag[0], np.vectorize(len)(frag[0])))
            for region in Region(frag=frag[0][sortIndex], raw=frag[1][sortIndex], depth=cov[sortIndex], \
                                 chr_=unit_chr_[sortIndex], start=unit_start[sortIndex], \
                                 end=unit_end[sortIndex], strand=unit_strand[sortIndex], \
                                 reads=unit_reads[sortIndex], cigars=unit_cigar[sortIndex], \
                                 unitStart=unit_rs[sortIndex], unitEnd=unit_re[sortIndex]).multi_fragments(
                unit_nfrag[sortIndex], frag[2][sortIndex]):
                result = SplitRegion(region).divide(single=False)
                for line in result.multi_printf():
                    loc.write(line)
                info, seq = result.writeInfo()
                out.write(f'{info}\n')
                fa.write(f'{seq}\n')
                for v in result.getVariants():
                    var.write(v)

def main():
    if not os.path.exists(args.output) and args.output != '.':
        os.makedirs(args.output, exist_ok=True)

    outUnit = f'{args.output}/unit.csv'
    outLoc = f'{args.output}/candidate_consolidated.csv'
    outInfo = f'{args.output}/final_eccDNA.csv'
    outSeq = f'{args.output}/consensus_sequence.fasta'
    outVar = f'{args.output}/variants.csv'
    with open(outUnit, 'w') as u:
        u.write('reads,Nfullpass,fragments,read_start,read_end,chr,start,end,strand,candidate_eccDNA,isGhost,order,cigar\n')
        for read in PAF(args.paf).reads():
            allcircle = read.bootstrap()
            if allcircle:
                for circle in allcircle:
                    if circle:
                        u.write(f'{circle.output_unit_info(read)}')
    printf(outUnit, outLoc, outInfo, outSeq, outVar)

if __name__ == '__main__':
    main()
