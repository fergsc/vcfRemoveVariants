#########################
## take a vcf as input
## randomly remove variants of a specified percent from 1, 2, ... all samples in vcf
## 1 sample has 10% removed, 2 samples have 10% removed, ..., 99 samples have 99% removed, 100 samples have 99% removed.
## Do this numRepeats times.
## written to test imputation accuracy.
#########################


library(vcfR)

vcfFileName = 'data/human_reduced.vcf'
removePC = c(0.90, 0.95, 0.97, 0.99)
numRepeats = 10

vcf = read.vcfR(vcfFileName, verbose = FALSE)
genos = vcf@gt

variants = dim(genos) 
samples = variants[2] - 1 # -1 as first column = "GT"
variants = variants[1]

write(fileName, "skip.txt", append=FALSE) # clear log file, stores the name of skipped/pre-existing vcf.

for(repeats in 1:numRepeats)
{
  for(missing in removePC)
  {
    numMissing = floor(variants*missing) # calcualte the number of variants to set to NA
    workingGenos = genos # reset working genotypes to default
    
    randSamples = sample(2:(samples+1), samples) # generate the list of samples to remove variants from, randomly
    
    currSample = 1
    for(s in randSamples)
    {
      fileName = sprintf("data/output/%i.%i_samples_%gpc.vcf.gz", repeats, currSample, missing*100)
      
      # incase this does not go to completion and needs to be restareted check for existing files.
      if(file.exists(fileName))
      {
        write(fileName, "skip.txt", append=TRUE)
        next
      }
      
      setNA = sample(1:variants, numMissing) # get list of random variants to set to NA
      workingGenos[setNA,s] = NA
      vcf@gt = workingGenos
      fileName = sprintf("data/output/%i.%i_samples_%gpc.vcf.gz", repeats, currSample, missing*100)
      print(fileName)
      write.vcf(vcf, fileName) 
      currSample = currSample+1
    }
  }
}
