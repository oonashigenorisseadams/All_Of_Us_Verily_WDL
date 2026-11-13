version 1.0

## aou_ld_panels_v9 — thinned LDstore2 matrices from the v5 panels
##
## Per group x chromosome:
##   prep:     v5 panel -> --mac 40 within the group -> BGEN 1.2 (8-bit, ref-first) -> bgenix index
##             -> variant table (index, ID, POS, REF, ALT, ALT_FREQ, OBS_CT)
##             -> windows: each 2 Mb "anchor" segment [a, b) gets the window [a, b + 2 Mb); anchors whose
##                window would hold more than max_window_variants are halved (down to min_anchor_bp),
##                because LDstore2 crashed on windows of ~64k+ variants. Anchors tile the chromosome,
##                so every pair within 2 Mb is still stored exactly once.
##             -> chunks of windows (each with its own BGEN slice)
##   ld_chunk: for each window, LDstore2 (--write-bcor --read-only-bgen --compression low, as in v3),
##             then keep pairs with r^2 > 0.005 and distance <= 2 Mb, storing each pair once
##             (left variant in the window's anchor); the dense .bcor is deleted right away
##   merge:    one Parquet of pairs (i, j, r) per group x chromosome
## r is signed, ALT-ALT positive, i/j index the variant table (0-based).
## Resuming: reuse_windows / reuse_pairs / reuse_stats (one inner array per group x chromosome, parallel to
## pgen_files) hand over finished 2 Mb-anchor windows from an earlier run; prep skips those anchors and
## merge adds their pair files, so only missing windows are computed.
## Requires window_bp = step_bp + max_dist_bp (default anchor 2 Mb + 2 Mb reach = 4 Mb windows).

workflow aou_ld_v5_matrices {
    input {
        # Parallel arrays: one entry per group x chromosome, e.g. .../anc_eur_75k.chr22.pgen
        Array[File] pgen_files
        Array[File] pvar_files
        Array[File] psam_files

        # optional, for resuming a partly finished run (see header)
        Array[Array[String]]? reuse_windows
        Array[Array[File]]?   reuse_pairs
        Array[Array[File]]?   reuse_stats

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
        Int    max_window_variants = 50000    # split anchors so LDstore2 never sees bigger windows
        Int    min_anchor_bp       = 250000

        String docker = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"

        Int prep_cpu             = 8
        Int prep_memory_gb       = 32
        Int prep_disk_multiplier = 25
        Int prep_disk_pad_gb     = 50

        Int ld_cpu          = 16
        Int ld_memory_gb    = 64     # peak was 75 GB at 66k variants; ~45 GB expected at 50k
        Int ld_disk_pad_gb  = 60
        Int ld_preemptible  = 2   # tries on spot VMs, then on-demand
    }

    scatter (i in range(length(pgen_files))) {
        Array[String] rw = if defined(reuse_windows) then select_first([reuse_windows])[i] else []
        Array[File]   rp = if defined(reuse_pairs)   then select_first([reuse_pairs])[i]   else []
        Array[File]   rs = if defined(reuse_stats)   then select_first([reuse_stats])[i]   else []

        call prep {
            input:
                pgen = pgen_files[i], pvar = pvar_files[i], psam = psam_files[i],
                plink2 = plink2, tools = tools,
                mac = mac, window_bp = window_bp, step_bp = step_bp, max_dist_bp = max_dist_bp,
                windows_per_chunk = windows_per_chunk,
                max_window_variants = max_window_variants, min_anchor_bp = min_anchor_bp,
                reuse_windows = rw,
                docker = docker, cpu = prep_cpu, memory_gb = prep_memory_gb,
                disk_multiplier = prep_disk_multiplier, disk_pad_gb = prep_disk_pad_gb
        }

        scatter (c in prep.chunks) {
            call ld_chunk {
                input:
                    chunk = c, ldstore = ldstore, tools = tools,
                    max_dist_bp = max_dist_bp, r2_min = r2_min,
                    docker = docker, cpu = ld_cpu, memory_gb = ld_memory_gb,
                    disk_pad_gb = ld_disk_pad_gb, preemptible = ld_preemptible
            }
        }

        call merge {
            input:
                parts = flatten(ld_chunk.pairs), thin_stats = ld_chunk.thin_stats, ld_stats = ld_chunk.ld_stats,
                reused_parts = rp, reused_stats = rs,
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
        Int max_window_variants
        Int min_anchor_bp
        Array[String] reuse_windows
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

        if [ $(( ~{window_bp} - ~{step_bp} )) -ne ~{max_dist_bp} ]; then
            echo "window_bp must equal step_bp + max_dist_bp" >&2; exit 1
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
        cat > make_windows.py <<'PYWIN'
#!/usr/bin/env python3
"""Adaptive LD windows for one chromosome.
Anchor segments [a, b) tile [0, last position]; each anchor gets the window [a, b + D), which holds
every pair whose left variant is in [a, b) and whose distance is <= D. An anchor is halved while its
window would hold more than nmax variants (down to min_anchor bp). Writes z/<window>.z (LDstore2 input)
for windows with >= 2 variants and windows.all.tsv:
name, start, end, first index, n variants, first id, anchor length."""
import os, sys
import numpy as np
var, label, S, D, nmax, amin = sys.argv[1], sys.argv[2], *map(int, sys.argv[3:7])
reuse = set(l.strip() for l in open(sys.argv[7]) if l.strip()) if len(sys.argv) > 7 else set()
rows = []
with open(var) as f:
    next(f)
    for line in f:
        t = line.rstrip("\n").split("\t")
        rows.append((t[1], t[2], int(t[3]), t[4], t[5]))      # id, chrom, pos, ref, alt
pos = np.fromiter((r[2] for r in rows), dtype=np.int64, count=len(rows))
if len(pos) == 0:
    sys.exit("empty variant table")
if np.any(np.diff(pos) < 0):
    sys.exit("variant table is not sorted by position")

def n_in(a, e):
    return int(np.searchsorted(pos, e, "left") - np.searchsorted(pos, a, "left"))

anchors = []
def add(a, b):
    if n_in(a, b + D) > nmax and (b - a) // 2 >= amin:
        m = a + (b - a) // 2
        add(a, m); add(m, b)
    else:
        anchors.append((a, b))
n_reused = 0
for k in range(int(pos[-1]) // S + 1):
    if f"{label}.{k * S}_{(k + 1) * S + D}" in reuse:   # finished in an earlier run
        n_reused += 1
        continue
    add(k * S, (k + 1) * S)

os.makedirs("z", exist_ok=True)
n_short = n_big = 0
with open("windows.all.tsv", "w") as out:
    for a, b in anchors:
        e = b + D
        lo, hi = int(np.searchsorted(pos, a, "left")), int(np.searchsorted(pos, e, "left"))
        n = hi - lo
        n_short += (b - a) < S
        n_big += n > nmax
        nm = f"{label}.{a}_{e}"
        out.write(f"{nm}\t{a}\t{e}\t{lo}\t{n}\t{rows[lo][0] if n else '.'}\t{b - a}\n")
        if n >= 2:
            with open(f"z/{nm}.z", "w") as z:
                z.write("rsid chromosome position allele1 allele2\n")
                for r in rows[lo:hi]:
                    z.write(f"{r[0]} {r[1]} {r[2]} {r[3]} {r[4]}\n")
print(f"reused from earlier run: {n_reused} (of {len(reuse)} given)  anchors to compute: {len(anchors)}  "
      f"shorter than {S} bp: {n_short}  still above {nmax} variants: {n_big}")
if n_reused != len(reuse):
    sys.exit("some reuse_windows names do not match this chromosome's 2 Mb anchors")
PYWIN
        tools/bin/python3 make_windows.py "${LABEL}.variants.tsv" "$LABEL" ~{step_bp} ~{max_dist_bp} \
            ~{max_window_variants} ~{min_anchor_bp} "~{write_lines(reuse_windows)}"

        # Drop windows with < 2 variants (no pairs).
        awk -F'\t' '$5 >= 2' windows.all.tsv > windows.tsv
        NWIN=$(wc -l < windows.tsv)
        echo "windows: $NWIN"

        echo "=== 5. chunks $(date -u +%FT%TZ) ==="
        # chunks = runs of <= windows_per_chunk adjacent windows (gaps left by reused windows start a new chunk)
        awk -F'\t' -v m=~{windows_per_chunk} '{ if (NR == 1 || k == m || $2 != nxt) { c++; k = 0; if (f) close(f); f = sprintf("chunklist.%04d", c - 1) }
                                               print > f; k++; nxt = $2 + $7 }' windows.tsv
        for CL in chunklist.*; do
            [ -e "$CL" ] || continue
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

        while IFS=$'\t' read -r NAME ST EN FIRST NV FID STEP; do
            echo "--- ${NAME} (${NV} variants) $(date -u +%FT%TZ)"
            printf 'z;bgen;bgi;sample;bdose;bcor;ld;n_samples\n' > master.txt
            printf '%s;%s;%s;%s;%s;%s;%s;%s\n' "work/z/${NAME}.z" work/chunk.bgen work/chunk.bgen.bgi work/panel.sample \
                w.bdose w.bcor w.ld "$N" >> master.txt
            for TRY in 1 2; do
                rm -f w.bcor w.bdose w.ld
                if "$TIME" -v -o w.time env LD_LIBRARY_PATH="$TL" ./ldstore --in-files master.txt --write-bcor --read-only-bgen \
                       --n-threads ~{cpu} --compression low > w.ldstore.log 2>&1; then break; fi
                echo "LDstore exited non-zero on ${NAME} (try ${TRY}):" >&2; tr '\r' '\n' < w.ldstore.log | tail -5 >&2
                if [ "$TRY" -eq 2 ]; then exit 1; fi
            done
            if [ ! -s w.bcor ] || grep -q '^Error' w.ldstore.log; then
                echo "LDstore did not produce a matrix for ${NAME}:" >&2; tr '\r' '\n' < w.ldstore.log | tail -20 >&2; exit 1
            fi
            WALL=$(awk 'index($0,"Elapsed (wall"){sub(/.*: /,""); print; exit}' w.time)
            RSS=$(awk 'index($0,"Maximum resident set size"){sub(/.*: /,""); printf "%.1f", $1/1048576; exit}' w.time)
            printf '%s\t%s\t%s\n' "$NAME" "$WALL" "$RSS" >> "${CH}.ld_stats.tsv"

            "$PY" thin_bcor.py --bcor w.bcor --start "$ST" --step "$STEP" --max-dist ~{max_dist_bp} \
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
        Array[File] reused_parts
        Array[File] reused_stats
        String label
        String tag
        File tools
        String docker
    }

    Int disk_gb = ceil((size(parts, "GB") + size(reused_parts, "GB")) * 2.5) + 20

    command <<<
        set -euo pipefail
        mkdir -p tools
        tar -xzf "~{tools}" -C tools
        cat > merge_pairs.py <<'PYMERGE'
#!/usr/bin/env python3
"""Concatenate per-window pair files (sorted by window start) into one Parquet per group x chromosome."""
import sys, re
import pyarrow.parquet as pq

out, listfile = sys.argv[1], sys.argv[2]
parts = [l.strip() for l in open(listfile) if l.strip()]
key = lambda p: int(re.search(r"\.(\d+)_\d+\.pairs\.parquet$", p).group(1))
parts = sorted(parts, key=key)
names = [re.sub(r".*/", "", p) for p in parts]
if len(set(names)) != len(names):
    sys.exit("duplicate window files")
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
        cat "~{write_lines(parts)}" "~{write_lines(reused_parts)}" > parts.txt
        echo "pair files: $(wc -l < parts.txt) (new $(wc -l < "~{write_lines(parts)}"), reused $(wc -l < "~{write_lines(reused_parts)}"))"
        tools/bin/python3 merge_pairs.py "~{label}.~{tag}.pairs.parquet" parts.txt

        # window, n_variants, pairs_checked, pairs_kept, thin_s, ldstore_wall, ldstore_max_rss_gb
        # (stats files from an earlier run can list windows that are not reused; keep only merged windows)
        export LC_ALL=C
        sed 's#.*/##; s#\.pairs\.parquet$##' parts.txt | sort -u > merged_windows.txt
        { cat "~{write_lines(thin_stats)}"; grep '\.thin_stats\.tsv$' "~{write_lines(reused_stats)}" || true; } | xargs -r cat \
            | sort -t$'\t' -k1,1 -u | join -t$'\t' - merged_windows.txt > thin.tsv
        { cat "~{write_lines(ld_stats)}"; grep '\.ld_stats\.tsv$' "~{write_lines(reused_stats)}" || true; } | xargs -r cat \
            | sort -t$'\t' -k1,1 -u > ld.tsv
        { printf 'window\tn_variants\tpairs_checked\tpairs_kept\tthin_s\tldstore_wall\tldstore_max_rss_gb\n'
          join -t$'\t' thin.tsv ld.tsv | sort -t. -k3,3n; } > "~{label}.~{tag}.window_stats.tsv"
        NS=$(( $(wc -l < "~{label}.~{tag}.window_stats.tsv") - 1 )); NM=$(wc -l < merged_windows.txt)
        echo "windows merged: $NM  with stats: $NS"
        if [ "$NS" -ne "$NM" ]; then echo "stats missing for some windows:" >&2; cut -f1 thin.tsv | comm -13 - merged_windows.txt | head >&2; fi
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
