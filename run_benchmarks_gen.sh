#!/bin/sh
# run_benchmarks_gen.sh — Run PROPACK generated-matrix scaling benchmarks
#                         and print Markdown result tables.
#
# Usage: ./run_benchmarks_gen.sh [BUILD_DIR]
#        BUILD_DIR defaults to "build".
#
# Prerequisites:
#   cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON
#   cmake --build build

set -e

BUILD="${1:-build}"

# Resolve to absolute path.
BUILD="$(cd "$BUILD" 2>/dev/null && pwd)" || {
    echo "error: build directory '${1:-build}' not found." >&2
    echo "Run:  cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON && cmake --build build" >&2
    exit 1
}

# Run one generated-matrix benchmark (no working directory needed).
#   run_one <prefix>
run_one() {
    _pfx=$1
    _exe="$BUILD/${_pfx}propack_benchmark_gen"
    [ -x "$_exe" ] || return 0
    "$_exe" 2>&1
}

# Collect output from all built precision variants.
collect() {
    run_one d
    run_one s
    run_one z
    run_one c
}

# ---- sanity check ----
_found=0
for _p in d s z c; do
    [ -x "$BUILD/${_p}propack_benchmark_gen" ] && _found=1
done
if [ "$_found" = 0 ]; then
    echo "error: no benchmark_gen executables found in $BUILD." >&2
    echo "Run:  cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON && cmake --build build" >&2
    exit 1
fi

echo "Running generated-matrix scaling benchmarks (this may take several minutes)..." >&2

# ---- run benchmarks and pipe into awk for table generation ----
collect | awk '
# ========== helpers ==========

function getval(line, key,    i, rest) {
    i = index(line, key "=")
    if (i == 0) return ""
    rest = substr(line, i + length(key) + 1)
    sub(/^ */, "", rest)
    sub(/ .*/, "", rest)
    return rest + 0
}

function fmt_wc(v) {
    if (v == 0) return "0"
    if (v < 0.0001) return sprintf("%.1e", v)
    if (v < 0.001)  return sprintf("%.5f", v)
    if (v < 0.01)   return sprintf("%.4f", v)
    return sprintf("%.3f", v)
}

function fmt_t(v) {
    if (v == 0) return "0"
    return sprintf("%.1e", v)
}

# ========== parsing ==========

# Track current precision prefix
/Scaling Benchmark/ {
    if (/Double/)    pfx = "d"
    if (/Single/)    pfx = "s"
    if (/Complex16/) pfx = "z"
    if (/Complex8/)  pfx = "c"
    seen[pfx] = 1
}

# Track current matrix configuration
/Generated random sparse/ {
    # Extract dimensions and nnz from the line
    match($0, /\(([0-9]+)x([0-9]+)/, _a)
    curM = getval($0, "")   # will parse below
    mattype = "sparse"
}

/Generated diagonal/ {
    mattype = "diag"
}

/Dimensions:/ {
    curM = getval($0, "m")
    curN = getval($0, "n")
    dimkey = curM "x" curN
}

/\*\*\* neig =/ { neig = $NF + 0 }

/Driver:/ { drv = $2 }

/Time\(s\)/ {
    k = pfx SUBSEP mattype SUBSEP dimkey SUBSEP drv SUBSEP neig
    wc_median[k] = getval($0, "median")
    # Track unique configs
    cfgkey = mattype SUBSEP dimkey SUBSEP neig
    if (!(cfgkey in cfg_seen)) {
        cfg_seen[cfgkey] = 1
        ncfg++
        cfg_mattype[ncfg] = mattype
        cfg_dim[ncfg] = dimkey
        cfg_m[ncfg] = curM
        cfg_n[ncfg] = curN
        cfg_neig[ncfg] = neig
    }
}

/nopx=/ {
    k = pfx SUBSEP mattype SUBSEP dimkey SUBSEP drv SUBSEP neig
    cnt[k, "nopx"]    = getval($0, "nopx")
    cnt[k, "nreorth"] = getval($0, "nreorth")
    cnt[k, "ndot"]    = getval($0, "ndot")
    cnt[k, "nbsvd"]   = getval($0, "nbsvd")
}

/nrestart=/ {
    k = pfx SUBSEP mattype SUBSEP dimkey SUBSEP drv SUBSEP neig
    cnt[k, "nrestart"] = getval($0, "nrestart")
    cnt[k, "nlandim"]  = getval($0, "nlandim")
    cnt[k, "nsing"]    = getval($0, "nsing")
}

/tmvopx=/ {
    k = pfx SUBSEP mattype SUBSEP dimkey SUBSEP drv SUBSEP neig
    tmr[k, "mvopx"]   = getval($0, "tmvopx")
    tmr[k, "getu0"]   = getval($0, "tgetu0")
    tmr[k, "lanbpro"] = getval($0, "tlanbpro")
    tmr[k, "bsvd"]    = getval($0, "tbsvd")
}

/tlansvd=/ {
    k = pfx SUBSEP mattype SUBSEP dimkey SUBSEP drv SUBSEP neig
    tmr[k, "lansvd"]  = getval($0, "tlansvd")
    tmr[k, "ritzvec"] = getval($0, "tritzvec")
    tmr[k, "restart"] = getval($0, "trestart")
    tmr[k, "reorth"]  = getval($0, "treorth")
}

# ========== table generation ==========

END {
    np = split("d s z c", P)

    lab["d"] = "double";    lab["s"] = "single"
    lab["z"] = "complex16"; lab["c"] = "complex8"

    for (i = 1; i <= np; i++) { p = P[i]
        bp[p] = p "lanbpro"
        sv[p] = p "lansvd"
        ir[p] = p "lansvd_irl"
    }

    # Count available precisions
    npresent = 0
    for (i = 1; i <= np; i++) { p = P[i]
        if (p in seen) npresent++
    }
    if (npresent == 0) exit

    # ---------- header ----------
    printf "## Generated-Matrix Scaling Benchmark Results\n\n"
    printf "**Precisions:**"
    sep = " "
    for (i = 1; i <= np; i++) { p = P[i]
        if (!(p in seen)) continue
        printf "%s%s", sep, lab[p]
        sep = ", "
    }
    printf ". Median of 3 reps. jobu/jobv = n (no vectors).\n\n"
    printf "---\n\n"

    # ======================================================
    #  Per-driver tables grouped by matrix type
    # ======================================================
    nd = split("lanbpro lansvd lansvd_irl", DRVBASE)

    nmt = split("sparse diag", MT)
    mtlab["sparse"] = "Random Sparse (nnz/col=20)"
    mtlab["diag"]   = "Diagonal (alpha=1.0)"

    for (mt_i = 1; mt_i <= nmt; mt_i++) { mtype = MT[mt_i]

        printf "### %s\n\n", mtlab[mtype]

        # -- wall-clock time table --
        printf "#### Wall-Clock Time (seconds, median)\n\n"

        # Header row
        printf "| m | n | neig |"
        for (i = 1; i <= np; i++) { p = P[i]
            if (!(p in seen)) continue
            printf " %s | %s | %s |", bp[p], sv[p], ir[p]
        }
        printf "\n|----:|----:|-----:|"
        for (i = 1; i <= np; i++) { p = P[i]
            if (!(p in seen)) continue
            printf "--------:|--------:|------------:|"
        }
        printf "\n"

        for (ci = 1; ci <= ncfg; ci++) {
            if (cfg_mattype[ci] != mtype) continue
            m = cfg_m[ci]; n = cfg_n[ci]; ne = cfg_neig[ci]
            dk = cfg_dim[ci]

            printf "| %d | %d | %d |", m, n, ne
            for (i = 1; i <= np; i++) { p = P[i]
                if (!(p in seen)) continue
                kb = p SUBSEP mtype SUBSEP dk SUBSEP bp[p] SUBSEP ne
                ks = p SUBSEP mtype SUBSEP dk SUBSEP sv[p] SUBSEP ne
                ki = p SUBSEP mtype SUBSEP dk SUBSEP ir[p] SUBSEP ne
                printf " %s |", (kb in wc_median) ? fmt_wc(wc_median[kb]) : "\342\200\224"
                printf " %s |", (ks in wc_median) ? fmt_wc(wc_median[ks]) : "\342\200\224"
                printf " %s |", (ki in wc_median) ? fmt_wc(wc_median[ki]) : "\342\200\224"
            }
            printf "\n"
        }
        printf "\n"

        # -- operation counts --
        printf "#### Operation Counts (last rep)\n\n"
        printf "| m | n | neig | Driver | nopx | nreorth | ndot | nbsvd | nrestart | nlandim |\n"
        printf "|----:|----:|-----:|--------|-----:|--------:|-----:|------:|---------:|--------:|\n"
        for (ci = 1; ci <= ncfg; ci++) {
            if (cfg_mattype[ci] != mtype) continue
            m = cfg_m[ci]; n = cfg_n[ci]; ne = cfg_neig[ci]
            dk = cfg_dim[ci]
            # Use double precision as representative (or first available)
            for (i = 1; i <= np; i++) { p = P[i]
                if (!(p in seen)) continue
                for (di = 1; di <= nd; di++) {
                    drvname = p DRVBASE[di]
                    kk = p SUBSEP mtype SUBSEP dk SUBSEP drvname SUBSEP ne
                    if (!((kk, "nopx") in cnt)) continue
                    printf "| %d | %d | %d | %s | %d | %d | %d | %d | %d | %d |\n",
                        m, n, ne, drvname,
                        cnt[kk,"nopx"], cnt[kk,"nreorth"],
                        cnt[kk,"ndot"], cnt[kk,"nbsvd"],
                        cnt[kk,"nrestart"], cnt[kk,"nlandim"]
                }
                break  # only first available precision for counts
            }
        }
        printf "\n"

        # -- timer breakdown --
        printf "#### Timer Breakdown (seconds, last rep)\n\n"
        printf "| m | n | neig | Driver | tmvopx | tgetu0 | tlanbpro | tbsvd | tlansvd | trestart | treorth |\n"
        printf "|----:|----:|-----:|--------|-------:|-------:|---------:|------:|--------:|---------:|--------:|\n"
        for (ci = 1; ci <= ncfg; ci++) {
            if (cfg_mattype[ci] != mtype) continue
            m = cfg_m[ci]; n = cfg_n[ci]; ne = cfg_neig[ci]
            dk = cfg_dim[ci]
            for (i = 1; i <= np; i++) { p = P[i]
                if (!(p in seen)) continue
                for (di = 1; di <= nd; di++) {
                    drvname = p DRVBASE[di]
                    kk = p SUBSEP mtype SUBSEP dk SUBSEP drvname SUBSEP ne
                    if (!((kk, "mvopx") in tmr)) continue
                    printf "| %d | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
                        m, n, ne, drvname,
                        fmt_t(tmr[kk,"mvopx"]), fmt_t(tmr[kk,"getu0"]),
                        fmt_t(tmr[kk,"lanbpro"]), fmt_t(tmr[kk,"bsvd"]),
                        fmt_t(tmr[kk,"lansvd"]), fmt_t(tmr[kk,"restart"]),
                        fmt_t(tmr[kk,"reorth"])
                }
                break  # only first available precision for timers
            }
        }
        printf "\n---\n\n"
    }
}
'
