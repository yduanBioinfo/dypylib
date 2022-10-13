# Annotation Objects

## Main contents

* Objects represents the annotations of a genome.
* Objects are independent of annotation file formats.
* Objects are defined at the basis of biology facts.
* The relationship between objects must been definded.
* Enable a flexible way to define new Objects.

## Development environment

Bioconda

## Data engine

Try gffutils at first. A relationship database might be used.

## Classes

### Base classes

* GENT: A genome element.

It should be a genome element without gaps. Usually, a GENT is defined by one line record. EXON/CDS are typical GENT.

* SGENT: A segment of genome element.

It should be a genome element with gaps. Usually, a SGENT is defined by multiple record lines. mRNA/Gene/Chromosome are typical SGENT. 

### Biological classes

* Genome
* Chromosome
* Gene
* Transcript
* Exon
* Cds

## Engine

* gffutils. Init Genome and other objects top to bottom

## Reference

1. [Description for GFF format](https://gmod.org/wiki/GFF2)