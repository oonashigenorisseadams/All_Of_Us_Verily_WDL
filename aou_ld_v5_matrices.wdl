version 1.0

## aou_ld_panels_v9 — thinned LDstore2 matrices from the v5 panels
##
## Per group x chromosome:
##   prep:     v5 panel -> --mac 40 within the group -> BGEN 1.2 (8-bit, ref-first) -> bgenix index
##             -> variant table (index, ID, POS, REF, ALT, ALT_FREQ, OBS_CT)
##             -> 4 Mb windows sliding every 2 Mb -> chunks of windows (each with its own BGEN slice)
##   ld_chunk: for each window, LDstore2 (--write-bcor --read-only-bgen --compression low, as in v3),
##             then keep pairs with r^2 > 0.005 and distance <= 2 Mb, storing each pair once
##             (left variant in the first half of the window); the dense .bcor is deleted right away
##   merge:    one Parquet of pairs (i, j, r) per group x chromosome
## r is signed, ALT-ALT positive, i/j index the variant table (0-based).
## Requires window_bp - step_bp >= max_dist_bp so every pair within max_dist is in some window.

workflow aou_ld_v5_matrices {
    input {
        # Parallel arrays: one entry per group x chromosome, e.g. .../anc_eur_75k.chr22.pgen
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files

        File plink2
        File ldstore        # ldstore_v2.0_x86_64
        File tools          # tar.gz of a conda env with bin/bgenix, bin/time, bin/python3 (+ numpy, pyarrow)

        String tag               = "mac40_r2gt0.005_d2Mb"
        Int    mac               = 40
        Int    window_bp         = 4000000
        Int    step_bp           = 2000000
        Int    max_dist_bp       = 2000000
        Float  r2_min            = 0.005
        Int    windows_per_chunk = 8

        String docker = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"

        Int prep_cpu             = 8
        Int prep_memory_gb       = 32
        Int prep_disk_multiplier = 25
        Int prep_disk_pad_gb     = 50

        Int ld_cpu          = 16
        Int ld_memory_gb    = 104
        Int ld_disk_pad_gb  = 60
        Int ld_preemptible  = 2   # tries on spot VMs, then on-demand
    }

    scatter (i in range(length(pgen_files))) {
        call prep {
            input:
                pgen = pgen_files[i], pvar = pvar_files[i], psam = psam_files[i],
                plink2 = plink2, tools = tools,
                mac = mac, window_bp = window_bp, step_bp = step_bp, max_dist_bp = max_dist_bp,
                windows_per_chunk = windows_per_chunk,
                docker = docker, cpu = prep_cpu, memory_gb = prep_memory_gb,
                disk_multiplier = prep_disk_multiplier, disk_pad_gb = prep_disk_pad_gb
        }

        scatter (c in prep.chunks) {
            call ld_chunk {
                input:
                    chunk = c, ldstore = ldstore, tools = tools,
                    step_bp = step_bp, max_dist_bp = max_dist_bp, r2_min = r2_min,
                    docker = docker, cpu = ld_cpu, memory_gb = ld_memory_gb,
                    disk_pad_gb = ld_disk_pad_gb, preemptible = ld_preemptible
            }
        }

        call merge {
            input:
                parts = flatten(ld_chunk.pairs), thin_stats = ld_chunk.thin_stats, ld_stats = ld_chunk.ld_stats,
                label = prep.label, tag = tag, tools = tools, docker = docker
        }
    }

    output {
        Array[File]        pairs        = merge.pairs
        Array[File]        window_stats = merge.window_stats
        Array[File]        variants     = prep.variants
        Array[Array[File]] prep_logs    = prep.logs
    }
}

task prep {
    input {
        File pgen
        File pvar
        File psam
        File plink2
        File tools
        Int mac
        Int window_bp
        Int step_bp
        Int max_dist_bp
        Int windows_per_chunk
        String docker
        Int cpu
        Int memory_gb
        Int disk_multiplier
        Int disk_pad_gb
    }

    Int disk_gb      = ceil(size(pgen, "GB") * disk_multiplier) + disk_pad_gb
    Int plink_mem_mb = memory_gb * 1024 - 4096

    command <<<
        set -euo pipefail
        echo "start: $(date -u +%FT%TZ)  nproc: $(nproc)"; df -h . | tail -1

        if [ $(( ~{window_bp} - ~{step_bp} )) -lt ~{max_dist_bp} ]; then
            echo "window_bp - step_bp must be >= max_dist_bp" >&2; exit 1
        fi
        if [ $(( ~{window_bp} % ~{step_bp} )) -ne 0 ]; then
            echo "window_bp must be a multiple of step_bp" >&2; exit 1
        fi

        mkdir -p tools logs bgen z chunks
        tar -xzf "~{tools}" -C tools
        TL="$PWD/tools/lib"; TIME="$PWD/tools/bin/time"
        cp "~{plink2}" ./plink2; chmod 750 ./plink2

        PGEN="~{pgen}"; PFX="${PGEN%.pgen}"; LABEL="$(basename "$PFX")"
        echo "$LABEL" > label.txt
        ls -lL "${PFX}.pgen" "${PFX}.pvar" "${PFX}.psam"

        echo "=== 1. MAC >= ~{mac} -> BGEN $(date -u +%FT%TZ) ==="
        "$TIME" -v -o "logs/${LABEL}.export.time" \
            ./plink2 --pfile "$PFX" --mac ~{mac} --export bgen-1.2 bits=8 ref-first --output-chr chrM \
                     --threads ~{cpu} --memory ~{plink_mem_mb} --out "bgen/${LABEL}"
        cp "bgen/${LABEL}.log" "logs/${LABEL}.export.plink.log"
        ./plink2 --pfile "$PFX" --mac ~{mac} --freq --output-chr chrM \
                 --threads ~{cpu} --memory ~{plink_mem_mb} --out "bgen/${LABEL}.mac"

        echo "=== 2. bgenix index + list $(date -u +%FT%TZ) ==="
        env LD_LIBRARY_PATH="$TL" tools/bin/bgenix -g "bgen/${LABEL}.bgen" -index
        env LD_LIBRARY_PATH="$TL" tools/bin/bgenix -g "bgen/${LABEL}.bgen" -list > "bgen/${LABEL}.list"

        echo "=== 3. variant table $(date -u +%FT%TZ) ==="
        # bgenix -list: alternate_ids rsid chromosome position number_of_alleles first_allele alternative_alleles
        awk 'BEGIN { k = 0 }
             FNR==NR { if (FNR==1) { for (i=1;i<=NF;i++) { if ($i=="ID") d=i; if ($i=="ALT_FREQS") f=i; if ($i=="OBS_CT") o=i }; next }
                       af[$d] = $f "\t" $o; next }
             /^#/ || $1 == "alternate_ids" { next }
             { if (!h) { print "index\tid\tchrom\tpos\tref\talt\talt_freq\tobs_ct"; h=1 }
               print (k+0) "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $7 "\t" af[$2]; k++ }' \
            "bgen/${LABEL}.mac.afreq" "bgen/${LABEL}.list" > "${LABEL}.variants.tsv"
        NVAR=$(( $(wc -l < "${LABEL}.variants.tsv") - 1 ))
        CHR=$(awk -F'\t' 'NR==2{print $3}' "${LABEL}.variants.tsv")
        echo "variants (MAC >= ~{mac}): $NVAR  chrom: $CHR"
        case "$CHR" in chr[0-9XYM]*) ;; *) echo "unexpected chromosome name '$CHR'" >&2; exit 1 ;; esac

        echo "=== 4. windows $(date -u +%FT%TZ) ==="
        # Each variant goes into every window [st, st+W) that contains it (st a multiple of S).
        # Windows are contiguous runs of the variant table, so first index + local index = global index.
        awk -v W=~{window_bp} -v S=~{step_bp} -v L="$LABEL" '
            BEGIN { idx = 0 }
            /^#/ || $1 == "alternate_ids" { next }
            {
                p = $4; k0 = int(p / S); nw = int(W / S)
                for (j = 0; j < nw; j++) {
                    st = (k0 - j) * S
                    if (st < 0) continue
                    nm = L "." st "_" (st + W)
                    f = "z/" nm ".z"
                    if (!(nm in n)) { print "rsid chromosome position allele1 allele2" > f; first[nm] = idx; fid[nm] = $2; start[nm] = st }
                    print $2, $3, $4, $6, $7 > f
                    n[nm]++
                }
                idx++
            }
            END { for (nm in n) printf "%s\t%d\t%d\t%d\t%d\t%s\n", nm, start[nm], start[nm] + W, first[nm], n[nm], fid[nm] }' \
            "bgen/${LABEL}.list" | sort -t$'\t' -k2,2n > windows.all.tsv

        # Drop windows with < 2 variants (no pairs).
        awk -F'\t' '$5 >= 2' windows.all.tsv > windows.tsv
        awk -F'\t' '$5 < 2 { print "z/" $1 ".z" }' windows.all.tsv | xargs -r rm -f
        NWIN=$(wc -l < windows.tsv)
        echo "windows: $NWIN"

        echo "=== 5. chunks $(date -u +%FT%TZ) ==="
        split -l ~{windows_per_chunk} -d -a 4 windows.tsv chunklist.
        for CL in chunklist.*; do
            C="${CL#chunklist.}"
            D="c_${C}"; mkdir -p "$D/z"
            FROM=$(awk -F'\t' 'NR==1{print ($2 < 1 ? 1 : $2)}' "$CL")
            TO=$(awk -F'\t' 'END{print $3 - 1}' "$CL")
            env LD_LIBRARY_PATH="$TL" tools/bin/bgenix -g "bgen/${LABEL}.bgen" -incl-range "${CHR}:${FROM}-${TO}" > "$D/chunk.bgen"
            env LD_LIBRARY_PATH="$TL" tools/bin/bgenix -g "$D/chunk.bgen" -index
            GOT=$(env LD_LIBRARY_PATH="$TL" tools/bin/bgenix -g "$D/chunk.bgen" -list 2>/dev/null | awk '!/^#/ && $0 !~ /alternate_ids/' | wc -l)
            WANT=$(awk -F'\t' -v a="$FROM" -v b="$TO" 'NR>1 && $4>=a && $4<=b' "${LABEL}.variants.tsv" | wc -l)
            echo "chunk ${C}: ${CHR}:${FROM}-${TO}  variants in slice ${GOT} (expected ${WANT})"
            if [ "$GOT" -ne "$WANT" ] || [ "$GOT" -lt 2 ]; then echo "slice ${C} incomplete" >&2; exit 1; fi
            cp "bgen/${LABEL}.sample" "$D/panel.sample"
            cp "$CL" "$D/windows.tsv"
            cut -f1 "$CL" | while read -r nm; do mv "z/${nm}.z" "$D/z/"; done
            tar -cf "chunks/${LABEL}.c${C}.tar" -C "$D" .
            rm -rf "$D"
        done
        echo "chunks: $(ls chunks | wc -l)"

        cp windows.tsv "logs/${LABEL}.windows.tsv"
        gzip -c "${LABEL}.variants.tsv" > "${LABEL}.variants.tsv.gz"
        echo "end: $(date -u +%FT%TZ)"; df -h . | tail -1
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: 0
    }

    output {
        String      label    = read_string("label.txt")
        File        variants = glob("*.variants.tsv.gz")[0]
        Array[File] chunks   = glob("chunks/*.tar")
        Array[File] logs     = glob("logs/*")
    }
}

task ld_chunk {
    input {
        File chunk
        File ldstore
        File tools
        Int step_bp
        Int max_dist_bp
        Float r2_min
        String docker
        Int cpu
        Int memory_gb
        Int disk_pad_gb
        Int preemptible
    }

    Int disk_gb = ceil(size(chunk, "GB") * 2) + disk_pad_gb

    command <<<
        set -euo pipefail
        echo "start: $(date -u +%FT%TZ)  nproc: $(nproc)"
        mkdir -p tools work pairs
        tar -xzf "~{tools}" -C tools
        TL="$PWD/tools/lib"; TIME="$PWD/tools/bin/time"; PY="$PWD/tools/bin/python3"
        cp "~{ldstore}" ./ldstore; chmod 750 ./ldstore
        tar -xf "~{chunk}" -C work
        CH="$(basename "~{chunk}" .tar)"

        cat > thin_bcor.py <<'PYTHIN'
#!/usr/bin/env python3
"""Thin one LDstore2 .bcor window to sparse pairs.

Keeps pairs (x < y) with r^2 > --r2, POS_y - POS_x <= --max-dist, and POS_x in
[--start, --start + --step). With windows of length W sliding by S and max-dist <= W - S,
that last rule stores every nearby pair exactly once across overlapping windows.
Writes a Parquet file with columns i, j (chromosome-level variant indices, matching the
variants table) and r (signed, ALT-ALT positive as written by LDstore2).
.bcor layout follows LDstore2 / PolyFun's reader (bcor1.1, strictly lower triangle).
"""
import argparse, struct, sys, time
import numpy as np

NBYTES = {0: 2, 1: 4, 2: 8, 3: 1}
DTYPE = {1: np.uint8, 2: np.uint16, 4: np.uint32, 8: np.uint64}
INT_NA = {1: 208, 2: 53248, 4: 3489660928, 8: 14987979559889010688}


def read_header(fn):
    with open(fn, "rb") as f:
        if f.read(7) != b"bcor1.1":
            sys.exit(f"{fn}: not an LDstore .bcor file")
        _fsize, = struct.unpack("<Q", f.read(8))
        n_samples, = struct.unpack("<I", f.read(4))
        n, = struct.unpack("<I", f.read(4))
        comp, = struct.unpack("<B", f.read(1))
        offset, = struct.unpack("<Q", f.read(8))
        ids, pos = [], np.empty(n, dtype=np.int64)
        for k in range(n):
            _lbuf, _idx = struct.unpack("<II", f.read(8))
            lr, = struct.unpack("<H", f.read(2)); ids.append(f.read(lr).decode())
            pos[k], = struct.unpack("<I", f.read(4))
            lc, = struct.unpack("<H", f.read(2)); f.read(lc)
            la, = struct.unpack("<I", f.read(4)); f.read(la)
            lb, = struct.unpack("<I", f.read(4)); f.read(lb)
    return n_samples, n, comp, offset, ids, pos


def thin(fn, start, step, max_dist, r2_min, first_index):
    n_samples, n, comp, offset, ids, pos = read_header(fn)
    nb = NBYTES[comp]
    na = DTYPE[nb](INT_NA[nb])
    scale = np.float64(2.0) ** -(8 * nb - 2)
    npairs = n * (n - 1) // 2
    mm = np.memmap(fn, dtype=DTYPE[nb], mode="r", offset=offset, shape=(npairs,)) if npairs else None

    lo = np.searchsorted(pos, start, "left")
    hi = np.searchsorted(pos, start + step, "left")
    ymax = np.searchsorted(pos, pos + max_dist, "right")
    I, J, R = [], [], []
    n_checked = 0
    for x in range(lo, hi):
        m = int(ymax[x]) - x - 1
        if m <= 0:
            continue
        base = npairs - (n - x) * (n - x - 1) // 2
        seg = np.asarray(mm[base:base + m])
        r = seg.astype(np.float64) * scale - 1.0
        keep = (seg != na) & (r * r > r2_min)
        js = np.nonzero(keep)[0]
        n_checked += m
        if js.size:
            I.append(np.full(js.size, first_index + x, dtype=np.uint32))
            J.append((first_index + x + 1 + js).astype(np.uint32))
            R.append(r[js].astype(np.float32))
    cat = lambda L, dt: np.concatenate(L) if L else np.empty(0, dtype=dt)
    return dict(n_samples=n_samples, n=n, ids=ids, pos=pos, n_checked=n_checked,
                i=cat(I, np.uint32), j=cat(J, np.uint32), r=cat(R, np.float32))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bcor", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--step", type=int, required=True)
    ap.add_argument("--max-dist", type=int, required=True)
    ap.add_argument("--r2", type=float, required=True)
    ap.add_argument("--first-index", type=int, required=True)
    ap.add_argument("--expect-first-id", default="")
    ap.add_argument("--expect-n", type=int, default=-1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--name", required=True)
    a = ap.parse_args()

    t0 = time.time()
    res = thin(a.bcor, a.start, a.step, a.max_dist, a.r2, a.first_index)
    if a.expect_n >= 0 and res["n"] != a.expect_n:
        sys.exit(f"{a.name}: bcor has {res['n']} variants, expected {a.expect_n}")
    if a.expect_first_id and res["n"] and res["ids"][0] != a.expect_first_id:
        sys.exit(f"{a.name}: first variant {res['ids'][0]} != expected {a.expect_first_id}")

    import pyarrow as pa, pyarrow.parquet as pq
    tbl = pa.table({"i": res["i"], "j": res["j"], "r": res["r"]})
    pq.write_table(tbl, a.out, compression="zstd")
    with open(a.stats, "a") as s:
        s.write(f"{a.name}\t{res['n']}\t{res['n_checked']}\t{len(res['r'])}\t{time.time() - t0:.1f}\n")


if __name__ == "__main__":
    main()
PYTHIN

        : > "${CH}.thin_stats.tsv"
        : > "${CH}.ld_stats.tsv"
        N=$(awk 'NR > 2' work/panel.sample | wc -l)

        while IFS=$'\t' read -r NAME ST EN FIRST NV FID; do
            echo "--- ${NAME} (${NV} variants) $(date -u +%FT%TZ)"
            printf 'z;bgen;bgi;sample;bdose;bcor;ld;n_samples\n' > master.txt
            printf '%s;%s;%s;%s;%s;%s;%s;%s\n' "work/z/${NAME}.z" work/chunk.bgen work/chunk.bgen.bgi work/panel.sample \
                w.bdose w.bcor w.ld "$N" >> master.txt
            "$TIME" -v -o w.time env LD_LIBRARY_PATH="$TL" ./ldstore --in-files master.txt --write-bcor --read-only-bgen \
                --n-threads ~{cpu} --compression low > w.ldstore.log 2>&1 || { tail -40 w.ldstore.log >&2; exit 1; }
            if [ ! -s w.bcor ] || grep -q '^Error' w.ldstore.log; then
                echo "LDstore did not produce a matrix for ${NAME}:" >&2; tr '\r' '\n' < w.ldstore.log | tail -20 >&2; exit 1
            fi
            WALL=$(awk 'index($0,"Elapsed (wall"){sub(/.*: /,""); print; exit}' w.time)
            RSS=$(awk 'index($0,"Maximum resident set size"){sub(/.*: /,""); printf "%.1f", $1/1048576; exit}' w.time)
            printf '%s\t%s\t%s\n' "$NAME" "$WALL" "$RSS" >> "${CH}.ld_stats.tsv"

            "$PY" thin_bcor.py --bcor w.bcor --start "$ST" --step ~{step_bp} --max-dist ~{max_dist_bp} \
                --r2 ~{r2_min} --first-index "$FIRST" --expect-first-id "$FID" --expect-n "$NV" \
                --out "pairs/${NAME}.pairs.parquet" --stats "${CH}.thin_stats.tsv" --name "$NAME"
            rm -f w.bcor w.bdose w.ld
        done < work/windows.tsv

        echo "=== ld stats (window, wall, max RSS GB) ==="; cat "${CH}.ld_stats.tsv"
        echo "=== thin stats (window, n, pairs checked, pairs kept, seconds) ==="; cat "${CH}.thin_stats.tsv"
        echo "end: $(date -u +%FT%TZ)"
    >>>

    runtime {
        docker: docker
        cpu: cpu
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: preemptible
    }

    output {
        Array[File] pairs      = glob("pairs/*.parquet")
        File        thin_stats = glob("*.thin_stats.tsv")[0]
        File        ld_stats   = glob("*.ld_stats.tsv")[0]
    }
}

task merge {
    input {
        Array[File] parts
        Array[File] thin_stats
        Array[File] ld_stats
        String label
        String tag
        File tools
        String docker
    }

    Int disk_gb = ceil(size(parts, "GB") * 2.5) + 20

    command <<<
        set -euo pipefail
        mkdir -p tools
        tar -xzf "~{tools}" -C tools
        cat > merge_pairs.py <<'PYMERGE'
#!/usr/bin/env python3
"""Concatenate per-window pair files (sorted by window start) into one Parquet per group x chromosome."""
import sys, re
import pyarrow.parquet as pq

out, parts = sys.argv[1], sys.argv[2:]
key = lambda p: int(re.search(r"\.(\d+)_\d+\.pairs\.parquet$", p).group(1))
parts = sorted(parts, key=key)
writer, n = None, 0
for p in parts:
    t = pq.read_table(p)
    if writer is None:
        writer = pq.ParquetWriter(out, t.schema, compression="zstd")
    if t.num_rows:
        writer.write_table(t)
        n += t.num_rows
if writer is None:
    sys.exit("no input parts")
writer.close()
print(f"{out}\t{len(parts)} windows\t{n} pairs")
PYMERGE
        tools/bin/python3 merge_pairs.py "~{label}.~{tag}.pairs.parquet" ~{sep=" " parts}

        # window, n_variants, pairs_checked, pairs_kept, thin_s, ldstore_wall, ldstore_max_rss_gb
        export LC_ALL=C
        cat ~{sep=" " thin_stats} | sort -t$'\t' -k1,1 > thin.tsv
        cat ~{sep=" " ld_stats} | sort -t$'\t' -k1,1 > ld.tsv
        { printf 'window\tn_variants\tpairs_checked\tpairs_kept\tthin_s\tldstore_wall\tldstore_max_rss_gb\n'
          join -t$'\t' thin.tsv ld.tsv | sort -t. -k3,3n; } > "~{label}.~{tag}.window_stats.tsv"
        awk -F'\t' 'NR>1{n++; c+=$3; k+=$4} END{printf "windows=%d pairs_checked=%d pairs_kept=%d (%.3f%%)\n", n, c, k, (c?100*k/c:0)}' \
            "~{label}.~{tag}.window_stats.tsv"
        ls -l "~{label}.~{tag}.pairs.parquet"
    >>>

    runtime {
        docker: docker
        cpu: 2
        memory: "8 GB"
        disks: "local-disk ~{disk_gb} SSD"
        preemptible: 0
    }

    output {
        File pairs        = "~{label}.~{tag}.pairs.parquet"
        File window_stats = "~{label}.~{tag}.window_stats.tsv"
    }
}

