version 1.0

## aou_ld_panels_v9 — panels v5: multiallelics split + indels
##
## Same design as v4 (scatter over chromosomes, loop panels inside each task). Changes from v4:
##   * multiallelic sites are split into one biallelic record per ALT allele (plink2 multiallelics=-)
##     instead of being dropped by --max-alleles 2
##   * indels are kept (snps_only default false)
##   * split alleles are trimmed of bases shared at the end of REF and ALT (AT>TT -> A>T; POS unchanged)
##   * every variant gets a unique ID chrN:POS:REF:ALT (--output-chr chrM keeps the 'chr' prefix); '*' alleles are removed; any remaining
##     duplicate IDs are removed (first kept) and counted
##   * filters (geno, hwe, mac) are applied AFTER splitting, i.e. to each allele separately
##   * mac defaults to 1 (drop alleles with no copies in the group); the export floor (MAC >= 40)
##     is applied later, in the matrix step, so these panels can serve other thresholds too
##
## Output names: anc_eur_75k.keep -> anc_eur_75k.chr22.{pgen,pvar,psam}

workflow ld_panel_filter_v5 {
    input {
        # Parallel arrays, same length and order. index i = one chromosome.
        Array[String] chroms
        Array[File]   pgen_files
        Array[File]   pvar_files
        Array[File]   psam_files

        # One keep file per panel (FID IID). Same keeps as v3/v4.
        Array[File]   keep_files

        File plink2

        # QC — computed within each panel, after splitting
        Boolean snps_only   = false
        Float   geno        = 0.05
        String  hwe         = "1e-50 keep-fewhet"
        Int     mac         = 1
        Float?  maf                     # leave empty: no frequency cutoff

        String docker          = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
        Int    cpu             = 8
        Int    memory_gb       = 32
        String disk_type       = "SSD"
        Int    disk_multiplier = 5      # temp split panel + final panel for one group at a time
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
        Float       geno
        String      hwe
        Int         mac
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

        PGEN_PATH="~{pgen}"
        PFX="${PGEN_PATH%.pgen}"
        ls -lL "${PFX}.pgen" "${PFX}.pvar" "${PFX}.psam"

        mkdir -p out out/freq out/lists tmp

        # Rewrites the split pvar: 5 columns, trimmed alleles, unique CHROM:POS:REF:ALT IDs,
        # and a list of '*' alleles to exclude. Genotype coding is unchanged.
        cat > trim_split_pvar.awk <<'AWK'
/^##/ { next }
/^#CHROM/ {
    for (i = 1; i <= NF; i++) { if ($i == "#CHROM") c = i; if ($i == "POS") p = i; if ($i == "REF") r = i; if ($i == "ALT") a = i }
    print "#CHROM\tPOS\tID\tREF\tALT"; next
}
{
    ref = $r; alt = $a
    while (length(ref) > 1 && length(alt) > 1 && substr(ref, length(ref)) == substr(alt, length(alt))) {
        ref = substr(ref, 1, length(ref) - 1); alt = substr(alt, 1, length(alt) - 1)
    }
    id = $c ":" $p ":" ref ":" alt
    if (alt == "*") print id > STAR
    print $c "\t" $p "\t" id "\t" ref "\t" alt
}
AWK

        MANIFEST="manifest.chr~{chrom}.tsv"
        printf 'chrom\tpanel\tn_keep\tn_samples_out\tn_split_alleles\tn_star\tn_dup_removed\tn_variants_out\tn_snp\tn_indel_mnp\tn_dup_ids\tmin_mac\tbytes_out\tseconds\n' > "$MANIFEST"

        while read -r K; do
            PANEL="$(basename "$K" .keep)"
            TMP="tmp/${PANEL}.chr~{chrom}"
            OUT="out/${PANEL}.chr~{chrom}"
            OUTF="out/freq/${PANEL}.chr~{chrom}"

            echo "=== panel ${PANEL} start $(date -u +%FT%TZ) ==="
            T0=$(date +%s)

            # Pass 1: this group's samples, every multiallelic site split into biallelic records.
            ./plink2 \
                --pfile "$PFX" \
                --keep "$K" \
                --make-pgen multiallelics=- \
                --output-chr chrM \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$TMP"

            # Trim alleles, set unique IDs, list '*' alleles.
            : > "${TMP}.star.txt"
            awk -v STAR="${TMP}.star.txt" -f trim_split_pvar.awk "${TMP}.pvar" > "${TMP}.trim.pvar"
            NSPLIT=$(awk '!/^#/' "${TMP}.trim.pvar" | wc -l)
            NSTAR=$(wc -l < "${TMP}.star.txt")
            EXCL=""
            if [ "$NSTAR" -gt 0 ]; then EXCL="--exclude ${TMP}.star.txt"; fi

            # Pass 2: per-allele QC within the group.
            ./plink2 \
                --pgen "${TMP}.pgen" \
                --pvar "${TMP}.trim.pvar" \
                --psam "${TMP}.psam" \
                $EXCL \
                --rm-dup force-first list \
                ~{true="--snps-only" false="" snps_only} \
                --max-alleles 2 \
                --geno ~{geno} \
                --hwe ~{hwe} \
                --mac ~{mac} \
                ~{"--maf " + maf} \
                --make-pgen \
                --output-chr chrM \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$OUT"

            NDUPRM=0
            if [ -f "${OUT}.rmdup.list" ]; then
                NDUPRM=$(wc -l < "${OUT}.rmdup.list")
                mv "${OUT}.rmdup.list" "out/lists/${PANEL}.chr~{chrom}.rmdup.list"
            fi
            cp "${TMP}.star.txt" "out/lists/${PANEL}.chr~{chrom}.star.txt"
            cp "${TMP}.log" "out/${PANEL}.chr~{chrom}.split.log"
            rm -f "${TMP}".*

            # Per-allele frequencies for the transfer package and for checks.
            ./plink2 \
                --pfile "$OUT" \
                --freq \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$OUTF"

            T1=$(date +%s)
            echo "=== panel ${PANEL} done  $(date -u +%FT%TZ) ==="

            NKEEP=$(awk '!/^#/ && NF' "$K" | wc -l)
            NSAMP=$(awk '!/^#/' "${OUT}.psam" | wc -l)
            NVAR=$(awk '!/^#/' "${OUT}.pvar" | wc -l)
            NSNP=$(awk '!/^#/ && length($4) == 1 && length($5) == 1' "${OUT}.pvar" | wc -l)
            NIND=$((NVAR - NSNP))
            NDUP=$(awk '!/^#/{print $3}' "${OUT}.pvar" | sort | uniq -d | wc -l)
            MINMAC=$(awk 'NR==1{for(i=1;i<=NF;i++){if($i=="ALT_FREQS")f=i; if($i=="OBS_CT")o=i}; next}
                          {a=$f*$o; m=($o-a<a)?($o-a):a; if(n==0||m<min)min=m; n++}
                          END{printf "%.0f", min}' "${OUTF}.afreq")
            BYTES=$(du -cb "${OUT}".pgen "${OUT}".pvar "${OUT}".psam | tail -1 | cut -f1)

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "~{chrom}" "$PANEL" "$NKEEP" "$NSAMP" "$NSPLIT" "$NSTAR" "$NDUPRM" "$NVAR" "$NSNP" "$NIND" "$NDUP" "$MINMAC" "$BYTES" "$((T1-T0))" \
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
        Array[File] out_lists = glob("out/lists/*")
        File        manifest  = "manifest.chr~{chrom}.tsv"
    }
}
