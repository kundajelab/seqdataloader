import argparse
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser(description="bin entire genome with pre-specified bin sizes and sliding window stride")
    parser.add_argument("--chrom_sizes",help="tab-delimited file with chromosome name in column 1 and chromosome size in column 2",default="hg19.chromsizes.txt")
    parser.add_argument("--bin_size",type=int,default=1000)
    parser.add_argument("--stride",type=int,default=50)
    parser.add_argument("--outf",default="hg19.bins.bed")
    return parser.parse_args()

def main():
    args=parse_args()
    chrom_sizes=pd.read_table(args.chrom_sizes,header=None,sep='\t')
    outf=open(args.outf,'w')

    for index,row in chrom_sizes.iterrows():
        chrom_name=row[0]
        chrom_size=row[1]
        print(chrom_name) 
        for bin_start in range(1,chrom_size-args.bin_size,args.stride):
            bin_end=bin_start+args.bin_size
            outf.write(chrom_name+'\t'+str(bin_start)+'\t'+str(bin_end)+'\n')
if __name__=="__main__":
    main()
    
