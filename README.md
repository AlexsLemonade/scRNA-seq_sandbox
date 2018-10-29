scRNA-seq_workflow

## Step 1: Identifying a dataset to use as an example pipeline:
### Q1: Which technology is most commonly used?
- 10xGenomics has been most common recently according to [Angerer et al, 2017](https://www.sciencedirect.com/science/article/pii/S245231001730077X)
- Smart-seq and it's paired-end sequencing version, smart-seq2 hae been freqently used over time. Many tumor datasets use one of these technologies or 10Xgenomics
- Parts of 10X's processing is properietary but otherwise it's processing is similar to smart-seq2
*Conclusion:* Smart-seq2 data would be good for this example
    
### Q2: Which tissues are most representative of medulloblastoma?
Medulloblastoma - subclass of primitive neuroectodermal tumor
    In order of relatedness (in my very unqualified opinion):
    1. Other neuroectodermal tumors
        - Ewing sarcomas (in the periphery, share a similar amount of neural differentiation)
        - Pineoblastomas (tumors of the pineal gland)
    2. Glial based tumors (but NOT microglia based) since these all branch off from neural cells at the same point
        - Gliomas
        - Oligodendrogliomas
        - Astrocytic gliomas
    
### Q3: What *datasets* are available at this point in time to represent these?
*Google search terms*: medulloblastoma/brain tumor single cell RNA-seq
*Top runners:*
- GSE70630: Tirosh et al, 2016 oligodendroglioma, n = 4825, smartseq2, FACS sorted
- GSE84465: Darmanis et al, 2017, glioblastoma, n = 3589, smart-seq2, plate sorted
- GSE102130: Tirosh et al, 2018, K27M mutant glioma, n = 4085, smart-seq2, plate sorted, RAW DATA NOT YET AVAILABLE
- GSE89567: Tirosh et al, 2017, IDH-mutant astrocytoma, n = 6341, smart-seq2, FACS sorting
- GSE57872: Anoop et al, 2014, primary glioblastoma, n = 875, smart-seq, plate sorted? 
 
 *Qualities to look for:*
 - Adequate sample size
 - Published data
 - Raw data available 
 - Data available in common format
 - Tissue relatedness to medulloblastoma 
 - Simple experimental design
 
## Step 2: Identifying a tools to use in pipeline
#### Questions for determining which tools to use:
- Which are applicable across many technologies? 
- Which are easy to use? 
- Which tools are well maintained and seem like they will be for the foreseeable future? 

*Preprocessing*: 
    sequence quality control: FASTQC
    adapter trimming: TrimGalore!, Prinseq
    genome alignment: STAR
    normalization: 
    
*Post-processing*: 
    Seurat
    ASAP
    Falco
    Scone
    scPipe (issue with this is that it is not applicable to paired end sequencing data)
