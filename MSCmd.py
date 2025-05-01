import click
import os
@click.group()
def cmd():
    pass

def removeExt(path):

    return os.path.splitext(path)[0]

from filter import FilterBed
from stats import StatsResult

@cmd.command()
@click.option('-bed', '--bed')
@click.option('-bed_out', '--bed_out')
@click.option('-source_path', '--source_path')
@click.option('-genome', '--genome')
def filter(bed, bed_out, source_path, genome="hg38"):

    FilterBed.filterBed(bed, bed_out, source_path, genome)
    statsfile = removeExt(bed_out)+"_stats.txt"
    StatsResult.stats_result(bed_out,statsfile)

from bam_manup import BamRecalib
@cmd.command()
@click.option('-bamin', '--bamin')
@click.option('-bamout', '--bamout')
@click.option('-filter_bed', '--filter_bed')
@click.option('-genome_ref', '--genome_ref')
def reflectToBam(bamin,bamout,filter_bed,genome_ref):

    BamRecalib.run_recalib(bamin,bamout,filter_bed,genome_ref)


from nnmodel import AttentionClassfication
@cmd.command()
@click.option('-base', '--base')
@click.option('-source_path', '--source_path')
@click.option('-genome', '--genome')
@click.option('-fp_ivtpath', '--fp_ivtpath')
@click.option('-outhistory', '--outhistory')
@click.option('-weightpath', '--weightpath')
def trainSequenceClassification(base,source_path, genome, fp_ivtpath,outhistory, weightpath):

    AttentionClassfication.trainNN(base,source_path, genome, fp_ivtpath, outhistory, weightpath, epoch=100)



if __name__ == '__main__': cmd()