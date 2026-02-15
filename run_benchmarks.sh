#!/bin/sh
# run_benchmarks.sh — Run PROPACK benchmarks and print Markdown result tables.
#
# Usage: ./run_benchmarks.sh [BUILD_DIR]
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

# Run one benchmark if the executable and data directory exist.
#   run_one <prefix> <type_dir>
run_one() {
    _pfx=$1 _tdir=$2
    _exe="$BUILD/${_pfx}propack_benchmark"
    _wdir="$BUILD/testing/$_tdir"
    [ -x "$_exe" ] && [ -d "$_wdir" ] || return 0
    ( cd "$_wdir" && "$_exe" ) 2>&1
}

# Collect output from all built precision variants.
collect() {
    run_one d double
    run_one s single
    run_one z complex16
    run_one c complex8
}

# ---- sanity check ----
_found=0
for _p in d s z c; do
    [ -x "$BUILD/${_p}propack_benchmark" ] && _found=1
done
if [ "$_found" = 0 ]; then
    echo "error: no benchmark executables found in $BUILD." >&2
    echo "Run:  cmake -B build -DPROPACK_BUILD_BENCHMARKS=ON && cmake --build build" >&2
    exit 1
fi

echo "Running benchmarks (this may take a minute)..." >&2

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

/PROPACK Double Precision/  { pfx = "d" }
/PROPACK Single Precision/  { pfx = "s" }
/PROPACK Complex16/         { pfx = "z" }
/PROPACK Complex8/          { pfx = "c" }

/Reading matrix from file/  { matfile[pfx] = $NF }

/Dimensions:/ {
    M[pfx] = getval($0, "m")
    N[pfx] = getval($0, "n")
}

/\*\*\* neig =/ { neig = $NF + 0 }

/Driver:/ { drv = $2 }

/Time\(s\)/ {
    k = pfx SUBSEP drv SUBSEP neig
    wc_median[k] = getval($0, "median")
}

/nopx=/ {
    k = pfx SUBSEP drv SUBSEP neig
    cnt[k, "nopx"]    = getval($0, "nopx")
    cnt[k, "nreorth"] = getval($0, "nreorth")
    cnt[k, "ndot"]    = getval($0, "ndot")
    cnt[k, "nbsvd"]   = getval($0, "nbsvd")
}

/nrestart=/ {
    k = pfx SUBSEP drv SUBSEP neig
    cnt[k, "nrestart"] = getval($0, "nrestart")
    cnt[k, "nlandim"]  = getval($0, "nlandim")
    cnt[k, "nsing"]    = getval($0, "nsing")
}

/tmvopx=/ {
    k = pfx SUBSEP drv SUBSEP neig
    tmr[k, "mvopx"]   = getval($0, "tmvopx")
    tmr[k, "getu0"]   = getval($0, "tgetu0")
    tmr[k, "lanbpro"] = getval($0, "tlanbpro")
    tmr[k, "bsvd"]    = getval($0, "tbsvd")
}

/tlansvd=/ {
    k = pfx SUBSEP drv SUBSEP neig
    tmr[k, "lansvd"]  = getval($0, "tlansvd")
    tmr[k, "ritzvec"] = getval($0, "tritzvec")
    tmr[k, "restart"] = getval($0, "trestart")
    tmr[k, "reorth"]  = getval($0, "treorth")
}

# ========== table generation ==========

END {
    np = split("d s z c", P)
    nn = split("10 50 100 200", NE)

    lab["d"] = "double";    lab["s"] = "single"
    lab["z"] = "complex16"; lab["c"] = "complex8"

    for (i = 1; i <= np; i++) { p = P[i]
        bp[p] = p "lanbpro"
        sv[p] = p "lansvd"
        ir[p] = p "lansvd_irl"
    }

    # ---------- header ----------
    printf "## Benchmark Results\n\n"
    printf "**Test matrices:**"
    sep = " "
    for (i = 1; i <= np; i++) { p = P[i]
        if (!(p in M)) continue
        printf "%s`%s` (%dx%d, %s)", sep, matfile[p], M[p], N[p], lab[p]
        sep = "; "
    }
    printf ". Median of 5 reps. jobu/jobv = n (no vectors).\n\n"
    printf "---\n\n"

    # ======================================================
    #  Xlanbpro — shared subroutine
    # ======================================================
    printf "### Lanczos Bidiagonalization (Xlanbpro) — shared subroutine\n\n"
    printf "`k = min(neig+10, kmax, min(m,n))`\n\n"

    # -- wall-clock time --
    printf "#### Wall-Clock Time (seconds, median)\n\n"
    printf "| neig |"
    for (i = 1; i <= np; i++) { p = P[i]
        if (p in M) printf " %s |", bp[p]
    }
    printf "\n|-----:|"
    for (i = 1; i <= np; i++) { p = P[i]
        if (p in M) printf "---------:|"
    }
    printf "\n"
    for (j = 1; j <= nn; j++) { ne = NE[j]
        printf "| %4d |", ne
        for (i = 1; i <= np; i++) { p = P[i]
            if (!(p in M)) continue
            k = p SUBSEP bp[p] SUBSEP ne
            printf " %s |", (k in wc_median) ? fmt_wc(wc_median[k]) : "\342\200\224"
        }
        printf "\n"
    }
    printf "\n"

    # -- operation counts --
    printf "#### Operation Counts (last rep)\n\n"
    printf "| neig | Prefix | nopx | nreorth | ndot |\n"
    printf "|-----:|--------|-----:|--------:|-----:|\n"
    for (j = 1; j <= nn; j++) { ne = NE[j]
        for (i = 1; i <= np; i++) { p = P[i]
            if (!(p in M)) continue
            k = p SUBSEP bp[p] SUBSEP ne
            if (!((k, "nopx") in cnt)) continue
            printf "| %4d | %s | %5d | %7d | %6d |\n",
                ne, p, cnt[k,"nopx"], cnt[k,"nreorth"], cnt[k,"ndot"]
        }
    }
    printf "\n"

    # -- timer breakdown --
    printf "#### Timer Breakdown (seconds, last rep)\n\n"
    printf "| neig | Prefix | tmvopx | tgetu0 | tlanbpro | treorth |\n"
    printf "|-----:|--------|-------:|-------:|---------:|--------:|\n"
    for (j = 1; j <= nn; j++) { ne = NE[j]
        for (i = 1; i <= np; i++) { p = P[i]
            if (!(p in M)) continue
            k = p SUBSEP bp[p] SUBSEP ne
            if (!((k, "mvopx") in tmr)) continue
            printf "| %4d | %s | %s | %s | %s | %s |\n",
                ne, p,
                fmt_t(tmr[k,"mvopx"]), fmt_t(tmr[k,"getu0"]),
                fmt_t(tmr[k,"lanbpro"]), fmt_t(tmr[k,"reorth"])
        }
    }
    printf "\n---\n\n"

    # ======================================================
    #  SVD Drivers
    # ======================================================
    printf "### SVD Drivers (Xlansvd, Xlansvd_irl)\n\n"
    printf "#### Wall-Clock Time (seconds, median)\n\n"

    # -- real precisions --
    nr = split("d s", RP)
    has_real = 0; for (i = 1; i <= nr; i++) if (RP[i] in M) has_real = 1

    if (has_real) {
        printf "| neig |"
        for (i = 1; i <= nr; i++) { p = RP[i]
            if (p in M) printf " %s | %s |", sv[p], ir[p]
        }
        printf "\n|-----:|"
        for (i = 1; i <= nr; i++) { p = RP[i]
            if (p in M) printf "--------:|------------:|"
        }
        printf "\n"
        for (j = 1; j <= nn; j++) { ne = NE[j]; printf "| %4d |", ne
            for (i = 1; i <= nr; i++) { p = RP[i]
                if (!(p in M)) continue
                ks = p SUBSEP sv[p] SUBSEP ne
                ki = p SUBSEP ir[p] SUBSEP ne
                printf " %s |", (ks in wc_median) ? fmt_wc(wc_median[ks]) : "\342\200\224"
                printf " %s |", (ki in wc_median) ? fmt_wc(wc_median[ki]) : "\342\200\224"
            }
            printf "\n"
        }
        printf "\n"
    }

    # -- complex precisions --
    nc = split("z c", CP)
    has_cpx = 0; for (i = 1; i <= nc; i++) if (CP[i] in M) has_cpx = 1

    if (has_cpx) {
        printf "| neig |"
        for (i = 1; i <= nc; i++) { p = CP[i]
            if (p in M) printf " %s | %s |", sv[p], ir[p]
        }
        printf "\n|-----:|"
        for (i = 1; i <= nc; i++) { p = CP[i]
            if (p in M) printf "--------:|------------:|"
        }
        printf "\n"
        for (j = 1; j <= nn; j++) { ne = NE[j]; printf "| %4d |", ne
            for (i = 1; i <= nc; i++) { p = CP[i]
                if (!(p in M)) continue
                ks = p SUBSEP sv[p] SUBSEP ne
                ki = p SUBSEP ir[p] SUBSEP ne
                printf " %s |", (ks in wc_median) ? fmt_wc(wc_median[ks]) : "\342\200\224"
                printf " %s |", (ki in wc_median) ? fmt_wc(wc_median[ki]) : "\342\200\224"
            }
            printf "\n"
        }
        printf "\n"
    }

    # -- operation counts --
    printf "#### Operation Counts (last rep)\n\n"
    printf "| neig | Driver | nopx | nreorth | ndot | nbsvd | nrestart | nlandim |\n"
    printf "|-----:|--------|-----:|--------:|-----:|------:|---------:|--------:|\n"
    for (i = 1; i <= np; i++) { p = P[i]
        if (!(p in M)) continue
        printf "| **%s (%dx%d)** | | | | | | | |\n", lab[p], M[p], N[p]
        for (j = 1; j <= nn; j++) { ne = NE[j]
            ks = p SUBSEP sv[p] SUBSEP ne
            ki = p SUBSEP ir[p] SUBSEP ne
            if ((ks, "nopx") in cnt)
                printf "| %4d | %s | %5d | %7d | %6d | %5d | %8d | %7d |\n",
                    ne, sv[p], cnt[ks,"nopx"], cnt[ks,"nreorth"],
                    cnt[ks,"ndot"], cnt[ks,"nbsvd"],
                    cnt[ks,"nrestart"], cnt[ks,"nlandim"]
            if ((ki, "nopx") in cnt)
                printf "| %4d | %s | %5d | %7d | %6d | %5d | %8d | %7d |\n",
                    ne, ir[p], cnt[ki,"nopx"], cnt[ki,"nreorth"],
                    cnt[ki,"ndot"], cnt[ki,"nbsvd"],
                    cnt[ki,"nrestart"], cnt[ki,"nlandim"]
        }
    }
    printf "\n"

    # -- timer breakdown --
    printf "#### Timer Breakdown (seconds, last rep)\n\n"
    printf "| neig | Driver | tmvopx | tgetu0 | tlanbpro | tbsvd | tlansvd | trestart | treorth |\n"
    printf "|-----:|--------|-------:|-------:|---------:|------:|--------:|---------:|--------:|\n"
    for (i = 1; i <= np; i++) { p = P[i]
        if (!(p in M)) continue
        printf "| **%s** | | | | | | | | |\n", lab[p]
        for (j = 1; j <= nn; j++) { ne = NE[j]
            ks = p SUBSEP sv[p] SUBSEP ne
            ki = p SUBSEP ir[p] SUBSEP ne
            if ((ks, "mvopx") in tmr)
                printf "| %4d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
                    ne, sv[p],
                    fmt_t(tmr[ks,"mvopx"]), fmt_t(tmr[ks,"getu0"]),
                    fmt_t(tmr[ks,"lanbpro"]), fmt_t(tmr[ks,"bsvd"]),
                    fmt_t(tmr[ks,"lansvd"]), fmt_t(tmr[ks,"restart"]),
                    fmt_t(tmr[ks,"reorth"])
            if ((ki, "mvopx") in tmr)
                printf "| %4d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
                    ne, ir[p],
                    fmt_t(tmr[ki,"mvopx"]), fmt_t(tmr[ki,"getu0"]),
                    fmt_t(tmr[ki,"lanbpro"]), fmt_t(tmr[ki,"bsvd"]),
                    fmt_t(tmr[ki,"lansvd"]), fmt_t(tmr[ki,"restart"]),
                    fmt_t(tmr[ki,"reorth"])
        }
    }
    printf "\n"
}
'
