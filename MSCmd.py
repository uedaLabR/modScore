def removeExt(path):
    return os.path.splitext(path)[0]

from filter import FilterBed
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
@click.option('-bamin', '--bamin')
@click.option('-bamout', '--bamout')
@click.option('-filter_bed', '--filter_bed')
def reflectToBam(bamin,bamout,filter_bed):

    BamRecalib.reflectToBam(bamin,bamout,filter_bed)



from nnmodel import TrainModel
@click.option('-fp_ivtpath', '--fp_ivtpath')
@click.option('-outhistory', '--outhistory')
@click.option('-source_path', '--source_path')
@click.option('-genome', '--genome')
def trainSequenceClassification(fp_ivtpath,outhistory,source_path, genome):


    TrainModel.trainNN(source_path, genome, fp_ivtpath, outhistory, eachsize=20000, epoch=100)