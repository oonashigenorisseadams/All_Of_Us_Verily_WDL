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
      --memory ~{mem * 900} --threads ~{threads} \
      --bfile ~{sub(bedfile, "\.bed$", "")} \
      --keep ~{pheno_file} \
      --mac 5 \
      --write-snplist \
      --out ~{outprefix}

    ~{plink2_file} \
      --memory ~{mem * 900} --threads ~{threads} \
      --bfile ~{sub(bedfile, "\.bed$", "")} \
      --keep ~{pheno_file} \
      --maf 0.1 \
      --pca approx \
      --out ~{outprefix}

    sed -i 's/^#//' ~{outprefix}.eigenvec

    # Join PCs onto the pheno/covar file by IID (column found by header name).
    # plink2 may omit the FID column from .eigenvec, so never join on FID+IID by position.
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR {
      if (FNR==1) {
        for (i=1; i<=NF; i++) if ($i=="IID") pi=i
        if (!pi) { print "ERROR: no IID column in pheno file header" > "/dev/stderr"; exit 1 }
        h1 = $0
        next
      }
      a[$pi] = $0
      next
    }
    FNR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="IID") ei=i
        else if ($i!="FID") pcs[++n]=i
      }
      if (!ei) { print "ERROR: no IID column in eigenvec header" > "/dev/stderr"; exit 1 }
      printf "%s", h1
      for (k=1; k<=n; k++) printf "%s%s", OFS, $pcs[k]
      printf "\n"
      next
    }
    ($ei in a) {
      printf "%s", a[$ei]
      for (k=1; k<=n; k++) printf "%s%s", OFS, $pcs[k]
      printf "\n"
    }
    ' ~{pheno_file} ~{outprefix}.eigenvec > ~{outprefix}_phenocovar.tsv

    # Fail loudly if the join matched nobody (header-only output)
    nrows=$(( $(wc -l < ~{outprefix}_phenocovar.tsv) - 1 ))
    echo "Samples in phenocovar after PC join: ${nrows}"
    if [ "${nrows}" -lt 1 ]; then
      echo "ERROR: PC join produced no samples; check that pheno IIDs match the .fam IIDs" >&2
      exit 1
    fi
  >>>

  output {
    File phenocovar = outprefix + "_phenocovar.tsv"
    File keeplist = outprefix + ".snplist"
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk 300 SSD"
  }
}

task step1 {
  input {
    File regenie_bin
    File phenocovar_file
    File keeplist
    File bedfile
    File bimfile
    File famfile
    String covariates
    String phenos
    String outprefix
    Int cpu
    Int mem
    Int threads
    String? catCovarList
  }

  command <<<
    set -euo pipefail

    echo "Staged files:"
    find . -maxdepth 6 -type f

    chmod 750 ~{regenie_bin}

    echo "regenie version check:"
    ~{regenie_bin} --version || true

    echo "Applying variant keeplist (--extract): ~{keeplist}"
    wc -l ~{keeplist}

    ~{regenie_bin} \
      --step 1 \
      --bed ~{sub(bedfile, "\.bed$", "")} \
      --ref-first \
      --extract ~{keeplist} \
      --phenoFile ~{phenocovar_file} \
      --phenoColList ~{phenos} \
      --covarFile ~{phenocovar_file} \
      --covarColList ~{covariates} \
      ~{"--catCovarList " + catCovarList} \
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
    docker: "gcc:12"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk 300 SSD"
  }
}

task step2 {
  input {
    File regenie_bin
    File bgen
    File sample_file
    File phenocovar_file
    Array[File] locos
    File pred
    String phenos
    String covariates
    Int min_mac
    Float min_info
    String outprefix_base
    String chr_name
    Int threads
    Int cpu
    Int mem
    String? catCovarList
  }
  
Int input_size_gb = ceil(size(bgen, "GB") + size(locos, "GB") + size(pred, "GB") + size(sample_file, "GB") + size(phenocovar_file, "GB"))
Int disk_size_gb  = input_size_gb * 2 + 50

  command <<<
    set -euo pipefail

    echo "Staged:"
    find . -maxdepth 6

    chmod 750 ~{regenie_bin}

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

    # pred.list may hold step 1's absolute paths, which don't exist in this task.
    # Rewrite each entry to point at the localized LOCO file in this working directory.
    while read -r pheno_name loco_path; do
      [ -z "${pheno_name}" ] && continue
      local_loco="$(pwd)/$(basename "${loco_path}")"
      if [ ! -e "${local_loco}" ]; then
        echo "ERROR: LOCO file for ${pheno_name} not localized: $(basename "${loco_path}")" >&2
        exit 1
      fi
      echo "${pheno_name} ${local_loco}"
    done < ~{pred} > pred_local.list

    echo "Rewritten pred list:"
    cat pred_local.list

    ~{regenie_bin} \
      --step 2 \
      --bgen ~{bgen} \
      --ref-first \
      --sample ~{sample_file} \
      --phenoFile ~{phenocovar_file} \
      --phenoColList ~{phenos} \
      --covarFile ~{phenocovar_file} \
      --covarColList ~{covariates} \
      ~{"--catCovarList " + catCovarList} \
      --pred pred_local.list \
      --bsize 400 \
      --minMAC ~{min_mac} \
      --minINFO ~{min_info} \
      --threads ~{threads} \
      --out ~{outprefix_base}_~{chr_name}
  >>>

  output {
    Array[File] assoc = glob(outprefix_base + "_" + chr_name + "_*.regenie")
    File log = outprefix_base + "_" + chr_name + ".log"
  }

  runtime {
    docker: "gcc:12"
    cpu: "~{cpu}"
    memory: "~{mem} GB"
    disks: "local-disk ~{disk_size_gb} SSD"
  }
}

workflow pca_then_regenie {
  input {
    File regenie_bin
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
    String? catCovarList
  }

  call pca_join_keeplist {
    input:
      pheno_file = pheno_file,
      plink2_file = plink2_file,
      bedfile = bedfile,
      bimfile = bimfile,
      famfile = famfile,
      outprefix = outprefix_pca,
      cpu = cpu_pca,
      mem = mem_pca,
      threads = threads_pca
  }

  call step1 {
    input:
      regenie_bin = regenie_bin,
      phenocovar_file = pca_join_keeplist.phenocovar,
      keeplist = pca_join_keeplist.keeplist,
      bedfile = bedfile,
      bimfile = bimfile,
      famfile = famfile,
      covariates = covariates,
      phenos = phenos,
      outprefix = outprefix_step1,
      cpu = cpu1,
      mem = mem1,
      threads = threads1,
      catCovarList = catCovarList
  }

  scatter (i in range(length(bgen_files))) {
    call step2 {
      input:
        regenie_bin = regenie_bin,
        bgen = bgen_files[i],
        sample_file = sample_files[i],
        phenocovar_file = pca_join_keeplist.phenocovar,
        locos = step1.locos,
        pred = step1.pred_file_raw,
        phenos = phenos,
        covariates = covariates,
        min_mac = min_mac,
        min_info = min_info,
        outprefix_base = outprefix_step2_base,
        chr_name = chr_names[i],
        threads = threads2,
        cpu = cpu2,
        mem = mem2,
        catCovarList = catCovarList
    }
  }
}
