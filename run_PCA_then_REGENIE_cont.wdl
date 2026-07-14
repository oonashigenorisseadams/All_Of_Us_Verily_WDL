version 1.0

task pca_join_keeplist {
  input {
    File plink2_file
    File pheno_file
    File bedfile
    File bimfile
    File famfile
    String outprefix
    Int cpu
    Int mem
    Int threads
  }

  command <<<
    set -euo pipefail
    echo "Staged files:"
    find . -maxdepth 6 -type f
    chmod 750 ~{plink2_file}
    ~{plink2_file} \
            --memory ~{mem}000 --threads ~{threads} \
            --bfile ~{sub(bedfile, "\.bed$", "")} \
            --keep ~{pheno_file} \
            --mac 5 \
            --write-snplist \
            --out ~{outprefix}
    ~{plink2_file} \
            --memory ~{mem}000 --threads ~{threads} \
            --bfile ~{sub(bedfile, "\.bed$", "")} \
            --keep ~{pheno_file} \
            --maf 0.1 \
            --pca approx \
            --out ~{outprefix}

    sed -i 's/^#//g' ~{outprefix}.eigenvec

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
    if (FNR==1) {
        h1 = $0
        next
    }
    key = $1 FS $2
    a[key] = $0
    next
}
FNR==1 {
    printf "%s", h1
    for (i=3; i<=NF; i++) printf "%s%s", OFS, $i
    printf "\n"
    next
}
{
    key = $1 FS $2
    if (key in a) {
        printf "%s", a[key]
        for (i=3; i<=NF; i++) printf "%s%s", OFS, $i
        printf "\n"
    }
}
' ~{pheno_file} ~{outprefix}.eigenvec > ~{outprefix}_phenocovar.tsv
  >>>

  output {
    File phenocovar = outprefix + "_phenocovar.tsv"
    File keeplist = outprefix + ".snplist"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk 300 HDD"
  }
}

task step1 {
  input {
    File phenocovar_file
    File bedfile
    File bimfile
    File famfile
    String covariates
    String phenos
    String outprefix
    Int cpu
    Int mem
    Int threads
  }

  command <<<
    set -euo pipefail
    echo "Staged files:"
    find . -maxdepth 6 -type f
        regenie \
            --step 1 \
            --bed ~{sub(bedfile, "\.bed$", "")} \
            --phenoFile ~{phenocovar_file} \
            --phenoColList ~{phenos} \
            --covarFile ~{phenocovar_file} \
            --covarColList ~{covariates} \
            --bsize 1000 \
            --lowmem \
            --threads ~{threads} \
            --out ~{outprefix}
  >>>

  output {
    File pred_file_raw = outprefix + "_pred.list"
    Array[File] locos = glob(outprefix + "_*.loco")
  }

  runtime {
    docker: "ghcr.io/erkkuleo/regenie-nonlinear:sha-ce562ea"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk 500 HDD"
    maxRetries: 3
  }
}

task step2 {
  input {
    File bgen
    File sample_file
    File phenocovar_file
    Array[File] locos
    File  pred
    String phenos
    String covariates
    Int min_mac
    Float min_info

    String outprefix_base
    String chr_name
    Int threads
    Int cpu
    Int mem
  }

  command <<<
  set -euo pipefail

  echo "Staged:"
  find . -maxdepth 6

  echo "Pred file contents:"
  cat ~{pred}

  echo "Localized LOCO files:"
  printf '%s\n' ~{sep=' ' locos}

  # Make the localized .loco files visible in the current working directory
  # using the exact basenames referenced inside the pred.list file.
  for f in ~{sep=' ' locos}; do
    bn=$(basename "$f")
    ln -sf "$f" "$bn"
  done

  echo "LOCO symlinks in cwd:"
  ls -lh *.loco || true

  regenie \
    --step 2 \
    --bgen ~{bgen} \
    --ref-first \
    --sample ~{sample_file} \
    --phenoFile ~{phenocovar_file} \
    --phenoColList ~{phenos} \
    --covarFile ~{phenocovar_file} \
    --covarColList ~{covariates} \
    --pred ~{pred} \
    --bsize 400 \
    --minMAC ~{min_mac} \
    --minINFO ~{min_info} \
    --threads ~{threads} \
    --out ~{outprefix_base}_~{chr_name}
  >>>

  output {
    Array[File] assoc = glob(outprefix_base + "_" + chr_name + "_*.regenie")
    File log   = outprefix_base + "_" + chr_name + ".log"
  }

  runtime {
    docker: "ghcr.io/erkkuleo/regenie-nonlinear:sha-ce562ea"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk 400 HDD"
  }
}

workflow pca_then_regenie {
  input {
    File pheno_file
    File plink2_file
    File bedfile
    File bimfile
    File famfile
    String covariates
    String phenos
    String outprefix_pca

    Int cpu_pca
    Int mem_pca
    Int threads_pca

    String outprefix_step1
    Int cpu1
    Int mem1
    Int threads1

    Array[File] bgen_files
    Array[File] sample_files
    Array[String] chr_names

    String outprefix_step2_base
    Int min_mac
    Float min_info
    Int cpu2
    Int mem2
    Int threads2
  }

  call pca_join_keeplist {
    input:
      pheno_file      = pheno_file,
      plink2_file     = plink2_file,
      bedfile         = bedfile,
      bimfile         = bimfile,
      famfile         = famfile,
      outprefix       = outprefix_pca,
      cpu             = cpu_pca,
      mem             = mem_pca,
      threads         = threads_pca
  }

  call step1 {
    input:
      phenocovar_file = pca_join_keeplist.phenocovar,
      bedfile         = bedfile,
      bimfile         = bimfile,
      famfile         = famfile,
      covariates      = covariates,
      phenos          = phenos,
      outprefix       = outprefix_step1,
      cpu             = cpu1,
      mem             = mem1,
      threads         = threads1
  }

  scatter (i in range(length(bgen_files))) {
    call step2 {
      input:
        bgen            = bgen_files[i],
        sample_file     = sample_files[i],
        phenocovar_file = pca_join_keeplist.phenocovar,
        locos           = step1.locos,
        pred            = step1.pred_file_raw,
        phenos          = phenos,
        covariates      = covariates,
        min_mac         = min_mac,
        min_info        = min_info,
        outprefix_base  = outprefix_step2_base,
        chr_name        = chr_names[i],
        threads         = threads2,
        cpu             = cpu2,
        mem             = mem2
    }
  }

}
