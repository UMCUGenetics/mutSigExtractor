# mutSigExtractor

An R packge for extracting SNV (COSMIC) signatures, indel signatures and SV signatures from vcf files.

The core functions for extracting signatures are:
```
## By default, returns the absolute contributions for the 30 COSMIC signatures
extractSigsSnv()

## Returns the counts of indels within repeat regions, indels with  flanking microhomology, and 
## indels not under the above 2 categories. Each category is further stratified by the length of 
## the indel.
extractSigsIndel()

## By default, returns the absolute contribution of the 6 SV signatures as described in this paper: 
## https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip
extractSigsSv()
```