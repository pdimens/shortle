genomesize = 900000000    # should be ~2x expected, for RAM purposes
readprefix = "Mme"
outprefix = "stonecrab"   # preferred output prefix
read_f = f"{readprefix}.F.fq.gz"
read_r = f"{readprefix}.R.fq.gz"
trim_f = f"{readprefix}.R1.fq.gz"
trim_r = f"{readprefix}.R2.fq.gz"

rule all:
    input:
        asm = f"polish_2/{outprefix}.hypo.fasta"
    threads: 20
    shell: "quast.py --large {input} --threads {threads} --conserved --min-contig 500 -o assembly_metrics"

rule trim:
    input:
        f = read_f,
        r = read_r
    output:
        f = trim_f,
        r = trim_r
    message: "Trimming with fastp"
    log:
        html = readprefix + "_trim_.html",
        txt = readprefix + "_trim_.txt"
    threads: 20
    params: "--cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 -l 50 -q 15 -u 50 --detect_adapter_for_pe"
    shell:
        """
        fastp --thread {threads} --in1 {input.f} --in2 {input.r} --out1 {output.f} --out2 {output.r} -j trim.json -h {log.html} {params} &> {log.txt}
        """

rule assemble_abyss:
    input:
        f = trim_f,
        r = trim_r
    output: expand(f"Abyss/{outprefix}" + "{i}",  i = ["-contigs.fa", "-unitigs.fa"])
    message: "Assembling using AbySS"
    log: expand(f"Abyss/{outprefix}-stats" + "{i}", i = [".csv", ".html", ".md", ".tab"])
    conda: "workflow/envs/assembly.yaml"
    threads: 20
    params:
        prefix = outprefix,
        gsize = f"G={genomesize}",
        other = "-C Abyss v=-v B=100G H=3 kc=5 j=10 k=81 "
    shell:
        """
    	mkdir -p tmp
        mkdir -p Abyss/tmp
        export TMPDIR=./tmp
        abyss-pe in='../{input.f} ../{input.r}' np={threads} name={params.prefix} {params.other} {params.gsize}
        """

rule redundans_1:
    input:
        contigs = f"Abyss/{outprefix}-contigs.fa",
        f = trim_f,
        r = trim_r
    output: "redundans_1/contigs.reduced.fa"
    conda: "workflow/envs/redundans.yaml"
    message: "Performing primary haplotig reduction and scaffolding with Redundans"
    threads: 20
    params: "-v -m 125 -o redundans_1 --usebwa --noscaffolding --nogapclosing --resume"
    shell:
        """
        redundans.py -i {input.f} {input.r} -f {input.contigs}  -t {threads} {params}
        """

rule decompress:
    input:
        f = trim_f,
        r = trim_r
    output:
        f = temp(outprefix + ".R1.fq"),
        r = temp(outprefix + ".R2.fq")
    message: "Decompressing short reads for SparseAssembler"
    shell:
        """
        gunzip -c {input.f} > {output.f}
        gunzip -c {input.r} > {output.r}
        """

rule sparseassembler:
    input:
        f = f"{outprefix}.R1.fq",
        r = f"{outprefix}.R2.fq"
    output: "SparseAssembler/Contigs.txt"
    message: "Assembling reads with SparseAssembler"
    params:
        other = "NodeCovTh 5 EdgeCovTh 3 k 90 g 15 RemoveChimera 1 ChimeraTh 1 ContigTh 1",
        gsize = "GS " + str(2*genomesize)

    shell:
        """
        mkdir -p SparseAssembler
        cd SparseAssembler
        SparseAssembler i1 ../{input.f} i2 ../{input.r} {params}
        cd ..
        """

rule dbg2olc:
    input:
        f = f"{outprefix}.R1.fq",
        r = f"{outprefix}.R2.fq",
        sparse = "SparseAssembler/Contigs.txt"
    output: "dbg2olc/DBG2OLC_Consensus.fasta"
    message: "Further Assembling Contigs from SparseAssembler using DBG2OLC"
    threads: 1
    params: "MinOverlap 115 PathCovTh 3 k 31 KmerCovTh 0"
    shell:
        """
        mkdir -p dbg2olc
        cd dbg2olc
        DBG2OLC Contigs ../{input.sparse} f ../{input.f} f ../{input.r} {params}
        cd ..
        """

rule map_reads_1:
    input:
        f = trim_f,
        r = trim_r,
        asm = "dbg2olc/DBG2OLC_Consensus.fasta"
    output: temp("mapped.bam")
    message: "Mapping reads to {input.asm}"
    threads: 20
    shell:
        """
		minimap2 -t {threads} -ax sr --secondary=no --MD --sam-hit-only {input.asm} {input.f} {input.r} | samtools view -hb -F4 -q5 -@ {threads} > {output}.tmp 2> /dev/null
		samtools sort -m 10G -l0  -@{threads} {output}.tmp > {output} && rm {output}.tmp
		samtools index {output} -@{threads}
        """

rule polish_1:
    input:
        f = trim_f,
        r = trim_r,
        bam = "mapped.bam",
        asm = "dbg2olc/DBG2OLC_Consensus.fasta"
    output: f"polish_1/{outprefix}.hypo.fasta"
    message: "Polishing {input.asm} with HyPo"
    threads: 10
    shell:
        """
        echo -e "{input.f}\n{input.r}" > shortread.files.txt
        hypo -d {input.asm} -r @shortread.files.txt -s 780m -c 100 -b {input.bam} -p 5 -t {threads} -o {output}
        """

rule map_reads2:
    input:
        f = trim_f,
        r = trim_r,
        asm = f"polish_1/{outprefix}.hypo.fasta"
    output:
        mapfile = temp("mapped2.bam")
    message: "Mapping reads to {input.asm}"
    threads: 20
    shell:
        """
        if [ ! -f "{output}.tmp" ]; then
	   	minimap2 -t {threads} -ax sr --secondary=no --MD --sam-hit-only {input.asm} {input.f} {input.r} | samtools view -hb -F4 -q5 -@ {threads} > {output}.tmp 2> /dev/null
        fi
		samtools sort -m 10G -l0  -@{threads} {output}.tmp > {output} && rm {output}.tmp
        samtools index {output} -@{threads}
        """

rule purge_haplotigs_histogram:
	input:
		asm = f"polish_1/{outprefix}.hypo.fasta",
		mapfile = "mapped2.bam"
	output:
		histo = "purge_haplotigs/mapped2.bam.gencov",
		histo_img = "purge_haplotigs/mapped2.bam.histogram.png"
	message: "Generating coverage histogram for purge_haplotigs"
	threads: 20
	params:
		depth = "-d 300"
	shell:
		"""
		mkdir -p purge_haplotigs
		cd purge_haplotigs/
		purge_haplotigs hist -b ../{input.mapfile} -g ../{input.asm} -t {threads} {params} 2> /dev/null
		"""

rule purge_haplotigs_suspects:
	input:
		hist_cov = "purge_haplotigs/mapped2.bam.gencov"
	output:
		cov_out = "purge_haplotigs/assembly_stats.csv"
	params:
		low = "-low 4",
		mid = "-mid 6",
		high = "-high 16"
	message: "Finding suspect contigs"
	shell:
		"""
		cd purge_haplotigs
		purge_haplotigs cov -i ../{input} {params} -o ../{output}
		"""


rule purge_haplotigs:
	input:
		asm = f"polish_1/{outprefix}.hypo.fasta",
		mapfile = "mapped2.bam",
		suspects = "purge_haplotigs/assembly_stats.csv"
	output:
		curated = f"purge_haplotigs/{outprefix}_purge.fasta"
	log:
		haplotigs = f"purge_haplotigs/{outprefix}_purge.haplotigs.fasta",
		artefacts = f"purge_haplotigs/{outprefix}_purge.artefacts.fasta",
		reassignments = f"purge_haplotigs/{outprefix}_purge.reassignments.tsv",
		logs = f"purge_haplotigs/{outprefix}_purge.contig_associations.log"
	threads: 20
	message: "Purging haplotigs"
	params:
		prefix = f"-o {outprefix}_purge"
	shell:
		"""
		cd purge_haplotigs
		purge_haplotigs purge {params} -t {threads} -g ../{input.asm} -c ../{input.suspects} -d -b ../{input.mapfile}
		"""


rule merge_assemblies:
    input:
        asm1 = "redundans_1/" + "contigs.reduced.fa",
        asm2 = "purge_haplotigs/" + outprefix + "_purge.fasta"
    output: f"merged/merged_{outprefix}.fasta"
    threads: 10
    params: f"-pre {outprefix}"
    shell: "cd merged; merge_wrapper.py {params} ../{input.asm1} ../{input.asm2}"

rule redundans_2:
    input:
        contigs = f"merged/merged_{outprefix}.fasta",
        f = trim_f,
        r = trim_r
    output: "redundans_2/" + "scaffolds.reduced.fa"
    conda: "workflow/envs/redundans.yaml"
    message: "Performing primary haplotig reduction on and scaffolding merged assembly with Redundans"
    threads: 20
    params: "-v -m 125 -o redundans_2 --usebwa --resume"
    shell:
        """
        redundans.py -i {input.f} {input.r} -f {input.contigs} -t {threads} {params}
        """

rule map_reads3:
    input:
        f = trim_f,
        r = trim_r,
        asm = "redundans_2/" + "scaffolds.reduced.fa"
    output: temp("polish_2/mapped3.bam")
    message: "Mapping reads to {input.asm}"
    threads: 20
    shell:
        """
        minimap2 --secondary=no --MD -ax sr -t {threads} {input.asm} {input.f} {input.r} | samtools view -Sb -F4 -q5 -@ {threads} - > {output}.tmp 2> /dev/null
    	samtools sort -m 16G -l0  -@{threads} {output}.tmp > {output} && rm {output}.tmp
        samtools index {output} -@{threads}
        """

use rule polish_1 as final_polish with:
    input:
        f = trim_f,
        r = trim_r,
        bam = "polish_2/" + "mapped3.bam",
        asm = "redundans_2/" + "scaffolds.reduced.fa"
    output: f"polish_2/{outprefix}.hypo.fasta"
