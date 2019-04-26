# mutSigExtractor

## Package description
mutSigExtractor is an R package for extracting SNV, indel and SV mutational signatures from vcf files. This performed in two main steps as will be described below.

### Counting mutation contexts
The first step involves counting the mutations belonging to specific contexts for each variant type:
- SNV: trinucleotide context, consisting of the point mutation and the 5' and 3' flanking nucleotides
- Indel: indels within repeat regions, indels with flanking microhomology; and other indels. Each category is further stratified by the repeat unit length, the number of bases in the indel sequence that are homologous, and the indel sequence length, respectively.
- SV: type (deletions, duplications, inversions, translocations) and length (0 to >10Mb)

### Determine signature contribution (by least squares fitting)
The contribution of each of the [30 COSMIC SNV signatures](https://cancer.sanger.ac.uk/cosmic/signatures) are then calculated from the SNV trinucleotide contexts using least squares fitting on the [signature profile matrix](https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt).

Similarly, the contribution of the SV signatures as described by [Nik-Zainal et al. 2016](https://www.nature.com/articles/nature17676) are calculated using the [SV signature profile matrix](https://media.nature.com/original/nature-assets/nature/journal/v534/n7605/extref/nature17676-s3.zip).

For indels, least squares fitting is performed. The contexts themselves serve as the mutational signatures.

## Getting started
The main functions for extracting signatures are:
```
extractSigsSnv()
extractSigsIndel()
extractSigsSv()
```

Note that SNVs and indels are often reported in the same vcf file. Therefore, extractSigsSnv() and extractSigsIndel() will automatically detect SNVs and indels, respectively (SNVs: REF length==1 and ALT length==1; indels: REF length>1 or ALT length>1). 

It is recommended that the vcf.filter argument is set to 'PASS' (or '.' for certain vcf files) to remove low quality variants. So for example:
```
extractSigsSnv('/path/to/vcf_with_snvs', vcf.filter='PASS')
```

While, by default, extractSigsSnv() and extractSigsSv() return mutational signature contributions, it is also possible to return the raw mutation context counts instead:
```
extractSigsSnv('/path/to/vcf_with_snvs', vcf.filter='PASS', output='contexts')
extractSigsSv('/path/to/vcf_with_svs', vcf.filter='PASS', output='contexts')
```

Ultimately, these functions return a one column data frame (essentially a vector) of the mutational signature contributions (or the mutation context counts), where the rownames are the names of the signatures/contexts, and the colname is the sample name (if provided). Using extractSigsIndel() as an example, the output will look like this:
```
extractSigsIndel('/path/to/vcf_with_indels', vcf.filter='PASS', sample.name='PD4115')

		PD4115
del.rep.len.1	224
del.rep.len.2	18
del.rep.len.3	17
del.rep.len.4	18
del.rep.len.5	17
ins.rep.len.1	152
ins.rep.len.2	26
ins.rep.len.3	4
ins.rep.len.4	2
ins.rep.len.5	29
del.mh.bimh.1	98
del.mh.bimh.2	185
del.mh.bimh.3	144
del.mh.bimh.4	84
del.mh.bimh.5	56
ins.mh.bimh.1	9
ins.mh.bimh.2	4
ins.mh.bimh.3	3
ins.mh.bimh.4	3
ins.mh.bimh.5	9
del.none.len.1	69
del.none.len.2	20
del.none.len.3	16
del.none.len.4	6
del.none.len.5	49
ins.none.len.1	18
ins.none.len.2	16
ins.none.len.3	5
ins.none.len.4	2
ins.none.len.5	10
```

With a small number of samples, it is possible to run these functions locally. With larger datasets however, it is advised to run these functions on an HPC. extractSigsIndel() is the most computationally demanding of the 3 functions, followed by extractSigsSnv(), and lastly extractSigsSv().


