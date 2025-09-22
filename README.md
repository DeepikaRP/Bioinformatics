<h2> Bacterial variant calling analysis </h2> 
<p> Variant calling is a DNA-Seq analysis to identify different variants and analyse like SNPs, InDels, SNVs, etc... Here we are trying to analyse the variants for E.coli bacteria treated under Doxycyline for AMR studies. </p>
<br>
<p> Workflow </p>
<p> QC : Fastqc </p>
<p> Reference genome De novo assembly : Spades or Velvet </p>
<p> Indexing & Alignment: BWA </p>
<p> Bam file indexing, sorting and filtering : Samtools</p>
<p> Variant calling : BCFTools</p>
<p> Annotation : SnpEff </p>
