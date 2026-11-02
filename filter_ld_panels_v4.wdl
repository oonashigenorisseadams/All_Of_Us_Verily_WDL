version 1.0

## aou_ld_panels_v9 — low-frequency panels (v4)
##
## Same design as filter_ld_panels.wdl (v3): scatter over chromosomes, loop panels
## inside each task. Differences from v3:
##   * frequency filter is configurable; default is count-based (--mac 20) instead of --maf 0.01
##   * HWE default relaxed to 1e-50 (keep-fewhet retained)
##   * --snps-only is a toggle (default on) so indels can be included without a code change
##   * writes per-panel allele frequencies (.afreq) and a min-MAC column in the manifest
##   * no merge step: LD matrices are computed per chromosome, so per-chromosome panels are
##     the deliverable
##
## Output names derive from basename(keep, ".keep"): anc_eur_75k.keep -> anc_eur_75k.chr22.{pgen,pvar,psam}

workflow ld_panel_filter_v4 {
    input {
        # Parallel arrays, same length and order. index i = one chromosome.
        Array[String] chroms
        Array[File]   pgen_files
        Array[File]   pvar_files
        Array[File]   psam_files

        # One keep file per panel (FID IID). Same keeps as v3.
        Array[File]   keep_files

        File plink2

        # QC — computed within each panel
        Boolean snps_only   = true
        Int     max_alleles = 2
        Float   geno        = 0.05
        String  hwe         = "1e-50 keep-fewhet"
        Int?    mac                     # e.g. 20; omit to disable
        Float?  maf                     # e.g. 0.01; omit to disable

        # Runtime. Pinned digest from job dd7cc235 (same as v3).
        String docker          = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
        Int    cpu             = 8
        Int    memory_gb       = 32
        String disk_type       = "SSD"
        Int    disk_multiplier = 4      # v3 used 3; outputs are larger without the MAF cutoff
        Int    disk_pad_gb     = 50
        Int    preemptible     = 0
    }

    scatter (i in range(length(chroms))) {
        call filter_chrom {
            input:
                chrom           = chroms[i],
                pgen            = pgen_files[i],
                pvar            = pvar_files[i],
                psam            = psam_files[i],
                keeps           = keep_files,
                plink2          = plink2,
                snps_only       = snps_only,
                max_alleles     = max_alleles,
                geno            = geno,
                hwe             = hwe,
                mac             = mac,
                maf             = maf,
                docker          = docker,
                cpu             = cpu,
                memory_gb       = memory_gb,
                disk_type       = disk_type,
                disk_multiplier = disk_multiplier,
                disk_pad_gb     = disk_pad_gb,
                preemptible     = preemptible
        }
    }

    output {
        Array[Array[File]] panel_pgen  = filter_chrom.out_pgen
        Array[Array[File]] panel_pvar  = filter_chrom.out_pvar
        Array[Array[File]] panel_psam  = filter_chrom.out_psam
        Array[Array[File]] panel_afreq = filter_chrom.out_afreq
        Array[Array[File]] plink_logs  = filter_chrom.out_logs
        Array[File]        manifests   = filter_chrom.manifest
    }
}

task filter_chrom {
    input {
        String      chrom
        File        pgen
        File        pvar
        File        psam
        Array[File] keeps
        File        plink2
        Boolean     snps_only
        Int         max_alleles
        Float       geno
        String      hwe
        Int?        mac
        Float?      maf
        String      docker
        Int         cpu
        Int         memory_gb
        String      disk_type
        Int         disk_multiplier
        Int         disk_pad_gb
        Int         preemptible
    }

    Int input_gb     = ceil(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))
    Int disk_gb      = input_gb * disk_multiplier + disk_pad_gb
    Int plink_mem_mb = memory_gb * 1024 - 4096

    command <<<
        set -euo pipefail

        echo "=== environment ==="
        echo "nproc:      $(nproc)"
        echo "mem:        $(awk '/MemTotal/{print $2/1048576" GB"}' /proc/meminfo)"
        echo "cwd device: $(df -h . | tail -1)"
        echo "start:      $(date -u +%FT%TZ)"

        cp "~{plink2}" ./plink2
        chmod 750 ./plink2
        ./plink2 --version

        # Prefix reconstruction. Batch mirrors the bucket layout, so the trio
        # co-locates and --pfile resolves all three.
        PGEN_PATH="~{pgen}"
        PFX="${PGEN_PATH%.pgen}"
        ls -lL "${PFX}.pgen" "${PFX}.pvar" "${PFX}.psam"

        mkdir -p out out/freq

        MANIFEST="manifest.chr~{chrom}.tsv"
        printf 'chrom\tpanel\tn_keep\tn_samples_out\tn_variants_out\tn_dup_ids\tmin_mac\tbytes_out\tseconds\n' > "$MANIFEST"

        # set -e: any panel failing aborts the task rather than leaving partial
        # outputs with rc=0.
        while read -r K; do
            PANEL="$(basename "$K" .keep)"
            OUT="out/${PANEL}.chr~{chrom}"
            OUTF="out/freq/${PANEL}.chr~{chrom}"

            echo "=== panel ${PANEL} start $(date -u +%FT%TZ) ==="
            T0=$(date +%s)

            ./plink2 \
                --pfile "$PFX" \
                --keep "$K" \
                ~{true="--snps-only" false="" snps_only} \
                --max-alleles ~{max_alleles} \
                --geno ~{geno} \
                --hwe ~{hwe} \
                ~{"--mac " + mac} \
                ~{"--maf " + maf} \
                --make-pgen \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$OUT"

            # Per-panel allele frequencies for the transfer package.
            ./plink2 \
                --pfile "$OUT" \
                --freq \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$OUTF"

            T1=$(date +%s)
            echo "=== panel ${PANEL} done  $(date -u +%FT%TZ) ==="

            # awk/wc rather than grep -c: grep exits 1 on zero matches and
            # would trip pipefail.
            NKEEP=$(awk '!/^#/ && NF' "$K" | wc -l)
            NSAMP=$(awk '!/^#/' "${OUT}.psam" | wc -l)
            NVAR=$(awk '!/^#/'  "${OUT}.pvar" | wc -l)
            NDUP=$(awk '!/^#/{print $3}' "${OUT}.pvar" | sort | uniq -d | wc -l)

            # Smallest minor allele count among retained variants (for the AoU
            # MAC >= 20 check). Columns located by header name.
            MINMAC=$(awk 'NR==1{for(i=1;i<=NF;i++){if($i=="ALT_FREQS")f=i; if($i=="OBS_CT")o=i}; next}
                          {a=$f*$o; m=($o-a<a)?($o-a):a; if(n==0||m<min)min=m; n++}
                          END{printf "%.0f", min}' "${OUTF}.afreq")

            BYTES=$(du -cb "${OUT}".pgen "${OUT}".pvar "${OUT}".psam | tail -1 | cut -f1)

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "~{chrom}" "$PANEL" "$NKEEP" "$NSAMP" "$NVAR" "$NDUP" "$MINMAC" "$BYTES" "$((T1-T0))" \
                >> "$MANIFEST"
        done < ~{write_lines(keeps)}

        echo "=== manifest ==="
        cat "$MANIFEST"
        echo "end: $(date -u +%FT%TZ)"
    >>>

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} ~{disk_type}"
        preemptible: preemptible
    }

    output {
        Array[File] out_pgen  = glob("out/*.pgen")
        Array[File] out_pvar  = glob("out/*.pvar")
        Array[File] out_psam  = glob("out/*.psam")
        Array[File] out_afreq = glob("out/freq/*.afreq")
        Array[File] out_logs  = glob("out/*.log")
        File        manifest  = "manifest.chr~{chrom}.tsv"
    }
}
