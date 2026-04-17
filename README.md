# Non-coding RNA Analysis

This repository contains scripts, data, and results for **non-coding RNA (ncRNA) analysis** as part of the M.Sc. Bioinformatics project.  

It includes preprocessing of RNA-seq data, genome indexing, alignment, transcript assembly, and functional annotation for ncRNAs.

---

## Project Overview
- **Goal:** Identify and characterize ncRNAs from RNA-Seq datasets.
- **Main Steps:**
  1. Quality control of raw reads (`fastp`)
  2. Genome index building (`hisat2`)
  3. Read alignment
  4. Transcript assembly and quantification
  5. Differential expression analysis
  6. Functional annotation (GO & KEGG enrichment)

---

## Repository Contents
- **Scripts/**
  - `extract_gene_info.awk` — pulls gene information
  - R scripts for differential expression & visualization
- **Data/**
  - Reference genomes, annotation files (`gencode`, `Homo_sapiens.GRCh38`)
  - Sample count matrices
- **Results/**
  - GO & KEGG enrichment results
  - Correlation analyses
- Large files such as `NONCODE_db.nsq` are managed via [Git LFS](https://git-lfs.github.com).

---

## Requirements
- **Tools:**
  - `fastp`
  - `hisat2`
  - `samtools`
  - `stringtie`
  - R (with `dplyr`, `ggplot2`, `clusterProfiler`)
- **Git LFS** for large files.

---

## Usage
Clone this repository (with LFS files):

