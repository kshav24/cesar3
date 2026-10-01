#!/usr/bin/env bash
# run_parity_suite.sh: Benchmark and verify exact output parity against baseline
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BASELINE_DIR="/home/keshav/cesar_baseline"

# The binary under test (CESAR 3.0).
BIN="${REPO_DIR}/cesar"

# The reference (CESAR 2.0) binary. Baselines MUST be produced by this, never by
# ${BIN} -- generating them with the binary under test only proves determinism,
# not parity, and would silently bless any regression.
BASELINE_BIN="${BASELINE_CESAR2:-/home/keshav/CESAR2.0/cesar}"

TEST_FILES=(
    "extra/example1.fa"
    "extra/example2.fa"
    "extra/example3.fa"
    "extra/example4.fa"
    "ENST00000262209.fa"
)

mkdir -p "${BASELINE_DIR}"

MODE="${1:-verify}" # "generate" or "verify"

if [ "${MODE}" = "generate" ]; then
    echo "=== GENERATING BASELINE OUTPUTS in ${BASELINE_DIR} ==="
    echo "=== reference binary: ${BASELINE_BIN} ==="
    if [ ! -x "${BASELINE_BIN}" ]; then
        echo "ERROR: reference (CESAR 2.0) binary not found or not executable: ${BASELINE_BIN}"
        echo "       Set BASELINE_CESAR2=/path/to/cesar and re-run."
        exit 1
    fi
    for tf in "${TEST_FILES[@]}"; do
        base_name="$(basename "${tf}")"
        out_file="${BASELINE_DIR}/${base_name}.out"
        echo -n "Running ${tf} -> ${out_file} ... "
        START_TIME=$(date +%s%N)
        "${BASELINE_BIN}" "${REPO_DIR}/${tf}" > "${out_file}"
        END_TIME=$(date +%s%N)
        ELAPSED=$(( (END_TIME - START_TIME) / 1000000 ))
        echo "DONE (${ELAPSED} ms, $(wc -l < "${out_file}") lines)"
    done
    echo "Baseline generation complete."
    exit 0
fi

echo "=== RUNNING PARITY AND TIMING SUITE against ${BASELINE_DIR} ==="
ALL_PASSED=1
TMP_OUT_DIR="$(mktemp -d /tmp/cesar_test_XXXXXX)"
trap 'rm -rf "${TMP_OUT_DIR}"' EXIT

for tf in "${TEST_FILES[@]}"; do
    base_name="$(basename "${tf}")"
    base_file="${BASELINE_DIR}/${base_name}.out"
    test_out="${TMP_OUT_DIR}/${base_name}.out"
    
    if [ ! -f "${base_file}" ]; then
        echo "ERROR: Baseline ${base_file} missing! Run with 'generate' first."
        exit 1
    fi
    
    START_TIME=$(date +%s%N)
    "${BIN}" "${REPO_DIR}/${tf}" > "${test_out}"
    END_TIME=$(date +%s%N)
    ELAPSED=$(( (END_TIME - START_TIME) / 1000000 ))
    
    if diff -u "${base_file}" "${test_out}" > "${TMP_OUT_DIR}/${base_name}.diff"; then
        echo "[PASS] ${tf}: 0 byte diff (parity 100%, ${ELAPSED} ms)"
    else
        echo "[FAIL] ${tf}: divergence detected! (${ELAPSED} ms)"
        head -n 20 "${TMP_OUT_DIR}/${base_name}.diff"
        ALL_PASSED=0
    fi
done

if [ ${ALL_PASSED} -eq 1 ]; then
    echo "=== ALL PARITY CHECKS PASSED ==="
    exit 0
else
    echo "=== PARITY VERIFICATION FAILED ==="
    exit 1
fi
