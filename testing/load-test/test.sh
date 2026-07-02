#!/bin/bash
set -euo pipefail

################################################################################
# WDROP Test Suite - Sequential Test Runner
#
# This script runs wdrop.sh multiple times sequentially with different modes and
# configurations to stress-test the pipeline.
#
# Usage: ./test.sh [OPTIONS]
#
# Options:
#   -u, --url URL              PDB URL to download
#                              (default: http://files.rcsb.org/download/1PSV.pdb)
#   -r, --reference-url URL    Optional reference PDB URL to download
#   -n, --num-tests NUM        Number of test iterations to run (default: 10)
#   -p, --path PATH            Base path for test directories (default: ./test_runs)
#   -m, --mode MODE            Mode to test: tinker, gromacs, or all (default: all)
#   -t, --type TYPE            Type to test: iteration count 1, 5, or all (default: all)
#   -l, --loop                 Run each mode/type combination sequentially in an
#                              infinite loop until stopped
#   -h, --help                 Show this help message
#
# Examples:
#   ./test.sh -u http://files.rcsb.org/download/1PSV.pdb -n 5 -p /tmp/tests -m all -t all
#   ./test.sh -u http://files.rcsb.org/download/1PSV.pdb -r http://files.rcsb.org/download/2PTC.pdb
#   ./test.sh -l -m gromacs -t 1
#
################################################################################

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
readonly SCRIPT="$PROJECT_ROOT/pipeline/wdrop.sh"

# Default values
PDB_URL="http://files.rcsb.org/download/1PSV.pdb"
REFERENCE_PDB_URL=""
NUM_TESTS=10
TEST_BASE_PATH="./test_runs"
TEST_MODES=("tinker" "gromacs")
TEST_TYPES=("1" "5")
LOOP_MODE=false

TEST_START_TIME=$(date +%s)

# Color codes for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

################################################################################
# Utility Functions
################################################################################

print_usage() {
  cat << EOF
Usage: $0 [OPTIONS]

Options:
  -u, --url URL              PDB URL to download
                             (default: http://files.rcsb.org/download/1PSV.pdb)
  -r, --reference-url URL    Optional reference PDB URL to download
  -n, --num-tests NUM        Number of test iterations to run
                             (default: 10)
  -p, --path PATH            Base path for test directories
                             (default: ./test_runs)
  -m, --mode MODE            Mode to test: tinker, gromacs, or all
                             (default: all)
  -t, --type TYPE            Type to test: iteration count 1, 5, or all
                             (default: all)
  -l, --loop                 Run each mode/type combination sequentially in an
                             infinite loop until stopped
  -h, --help                 Show this help message

Examples:
  # Default run with 10 tests
  ./test.sh

  # Run 5 tests with gromacs mode only
  ./test.sh -n 5 -m gromacs

  # Run single-cycle (1 iteration) tests only in /tmp/test_dir
  ./test.sh -p /tmp/test_dir -t 1

  # Run with an optional reference PDB
  ./test.sh -u http://files.rcsb.org/download/1PSV.pdb \
            -r http://files.rcsb.org/download/2PTC.pdb

  # Run sequential tests forever until interrupted
  ./test.sh -l -m all -t all

Notes:
  Standard mode:
    With -n 2 -t 1:
      - 2 tests × 2 modes (tinker, gromacs) = 4 total test runs
      - Tests run one after another
      - Each will create a separate folder with its own PDB download
      - If --reference-url is set, the reference PDB is also downloaded per test

  Loop mode (-l):
      - Runs tests one after another
      - Repeats forever until stopped with Ctrl+C or SIGTERM
      - Prints the runtime of each individual run
EOF
  exit 0
}

log_msg() {
  local level="$1"
  shift
  local msg="$*"
  local timestamp
  timestamp=$(date '+%Y-%m-%d %H:%M:%S')

  case "$level" in
    INFO)
      echo -e "${BLUE}[${timestamp}]${NC} ${msg}"
      ;;
    SUCCESS)
      echo -e "${GREEN}[${timestamp}]${NC} ✓ ${msg}"
      ;;
    ERROR)
      echo -e "${RED}[${timestamp}]${NC} ✗ ${msg}" >&2
      ;;
    WARN)
      echo -e "${YELLOW}[${timestamp}]${NC} ⚠ ${msg}"
      ;;
    *)
      echo "[${timestamp}] ${msg}"
      ;;
  esac
}

format_duration() {
  local total_seconds="$1"
  local hours=$((total_seconds / 3600))
  local minutes=$(((total_seconds % 3600) / 60))
  local seconds=$((total_seconds % 60))

  if (( hours > 0 )); then
    printf '%dh %dm %ds' "$hours" "$minutes" "$seconds"
  elif (( minutes > 0 )); then
    printf '%dm %ds' "$minutes" "$seconds"
  else
    printf '%ds' "$seconds"
  fi
}

cleanup() {
  log_msg WARN "Received interrupt. Stopping test runner..."
  exit 130
}

check_prerequisites() {
  log_msg INFO "Checking prerequisites..."

  if [[ ! -f "$SCRIPT" ]]; then
    log_msg ERROR "wdrop.sh not found at: $SCRIPT"
    exit 1
  fi

  if ! command -v curl &> /dev/null && ! command -v wget &> /dev/null; then
    log_msg ERROR "Neither curl nor wget found. Install one to download PDB files."
    exit 1
  fi

  log_msg SUCCESS "Prerequisites check passed"
}

parse_arguments() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      -u|--url)
        PDB_URL="$2"
        shift 2
        ;;
      -r|--reference-url)
        REFERENCE_PDB_URL="$2"
        shift 2
        ;;
      -n|--num-tests)
        NUM_TESTS="$2"
        shift 2
        ;;
      -p|--path)
        TEST_BASE_PATH="$2"
        shift 2
        ;;
      -m|--mode)
        if [[ "$2" == "all" ]]; then
          TEST_MODES=("tinker" "gromacs")
        else
          TEST_MODES=("$2")
        fi
        shift 2
        ;;
      -t|--type)
        if [[ "$2" == "all" ]]; then
          TEST_TYPES=("1" "5")
        else
          TEST_TYPES=("$2")
        fi
        shift 2
        ;;
      -l|--loop)
        LOOP_MODE=true
        shift
        ;;
      -h|--help)
        print_usage
        ;;
      *)
        log_msg ERROR "Unknown option: $1"
        print_usage
        ;;
    esac
  done
}

download_pdb() {
  local url="$1"
  local output="$2"

  if command -v curl &> /dev/null; then
    curl -f -s -L "$url" -o "$output" 2>/dev/null || return 1
  elif command -v wget &> /dev/null; then
    wget -q "$url" -O "$output" 2>/dev/null || return 1
  else
    return 1
  fi

  [[ -f "$output" ]]
}

extract_pdb_id() {
  local url="$1"
  basename "$url" .pdb
}

create_test_directory() {
  local test_num="$1"
  local mode="$2"
  local test_type="$3"
  local suffix=""

  if [[ "$LOOP_MODE" == true ]]; then
    suffix="_$(date +%Y%m%d_%H%M%S)"
  fi

  local test_dir="${TEST_BASE_PATH}/test_${test_num}_${mode}_${test_type}${suffix}"

  mkdir -p "$test_dir"
  echo "$test_dir"
}

run_single_test_body() {
  local test_num="$1"
  local mode="$2"
  local test_type="$3"
  local pdb_id
  pdb_id=$(extract_pdb_id "$PDB_URL")

  local reference_pdb_id=""
  if [[ -n "$REFERENCE_PDB_URL" ]]; then
    reference_pdb_id=$(extract_pdb_id "$REFERENCE_PDB_URL")
  fi

  local test_dir
  test_dir=$(create_test_directory "$test_num" "$mode" "$test_type")
  local pdb_file="${test_dir}/${pdb_id}.pdb"
  local reference_pdb_file=""
  if [[ -n "$REFERENCE_PDB_URL" ]]; then
    reference_pdb_file="${test_dir}/${reference_pdb_id}.pdb"
  fi

  local test_log="${test_dir}/test.log"
  local result_file="${test_dir}/result.txt"

  log_msg INFO "Starting test #${test_num}: MODE=$mode, TYPE=$test_type in $test_dir"

  if ! download_pdb "$PDB_URL" "$pdb_file"; then
    log_msg ERROR "Test #${test_num}: Failed to download PDB"
    echo "FAILED: Download error for input PDB" > "$result_file"
    return 1
  fi

  if [[ -n "$REFERENCE_PDB_URL" ]]; then
    if ! download_pdb "$REFERENCE_PDB_URL" "$reference_pdb_file"; then
      log_msg ERROR "Test #${test_num}: Failed to download reference PDB"
      echo "FAILED: Download error for reference PDB" > "$result_file"
      return 1
    fi
  fi

  log_msg INFO "Test #${test_num}: PDB(s) downloaded, starting wdrop.sh"

  local test_start
  test_start=$(date +%s)

  local cmd=("$SCRIPT" "$pdb_file" "$mode" "--iterations" "$test_type")
  if [[ -n "$REFERENCE_PDB_URL" ]]; then
    cmd+=("$reference_pdb_file")
  fi

  if "${cmd[@]}" >> "$test_log" 2>&1; then
    local test_end
    test_end=$(date +%s)
    local test_duration=$((test_end - test_start))
    local formatted_duration
    formatted_duration=$(format_duration "$test_duration")

    log_msg SUCCESS "Test #${test_num} ($mode/$test_type): Completed in ${formatted_duration}"
    echo "SUCCESS: Completed in ${formatted_duration}" > "$result_file"
    return 0
  else
    local exit_code=$?
    local test_end
    test_end=$(date +%s)
    local test_duration=$((test_end - test_start))
    local formatted_duration
    formatted_duration=$(format_duration "$test_duration")

    log_msg ERROR "Test #${test_num} ($mode/$test_type): Failed with exit code $exit_code after ${formatted_duration}"
    echo "FAILED: Exit code $exit_code after ${formatted_duration}" > "$result_file"

    if [[ -f "$test_log" ]]; then
      log_msg WARN "Last 20 lines of ${test_log}:"
      tail -20 "$test_log"
    fi

    return 1
  fi
}

run_single_test_sequential() {
  local test_num="$1"
  local mode="$2"
  local test_type="$3"
  local run_start
  run_start=$(date +%s)

  if {
    run_single_test_body "$test_num" "$mode" "$test_type"
  } 2>&1 | tee -a "${TEST_BASE_PATH}/global.log"; then
    local run_end
    run_end=$(date +%s)
    local run_duration=$((run_end - run_start))
    log_msg INFO "Run runtime: $(format_duration "$run_duration")"
    return 0
  else
    local status=$?
    local run_end
    run_end=$(date +%s)
    local run_duration=$((run_end - run_start))
    log_msg WARN "Run runtime: $(format_duration "$run_duration")"
    return "$status"
  fi
}

collect_results() {
  log_msg INFO "Collecting test results..."
  echo ""

  local total=0
  local passed=0
  local failed=0
  declare -a failed_tests=()

  for test_dir in "${TEST_BASE_PATH}"/test_*/; do
    [[ -d "$test_dir" ]] || continue

    ((total++))
    local test_name
    test_name=$(basename "$test_dir")
    local result_file="${test_dir}/result.txt"

    if [[ -f "$result_file" ]]; then
      local result
      result=$(cat "$result_file")

      if [[ "$result" =~ SUCCESS ]]; then
        ((passed++))
        log_msg SUCCESS "$test_name: $result"
      else
        ((failed++))
        log_msg ERROR "$test_name: $result"
        failed_tests+=("$test_name: $result")
      fi
    else
      ((failed++))
      log_msg ERROR "$test_name: No result file"
      failed_tests+=("$test_name: No result file")
    fi
  done

  echo ""
  log_msg INFO "======================================"
  log_msg INFO "Test Summary"
  log_msg INFO "======================================"
  log_msg INFO "Total tests: $total"
  log_msg SUCCESS "Passed: $passed"
  log_msg ERROR "Failed: $failed"
  log_msg INFO "======================================"

  if [[ ${#failed_tests[@]} -gt 0 ]]; then
    echo ""
    log_msg WARN "Failed tests details:"
    for fail_info in "${failed_tests[@]}"; do
      log_msg ERROR "  $fail_info"
    done
  fi

  echo ""
  return $failed
}

print_test_configuration() {
  log_msg INFO "======================================"
  log_msg INFO "WDROP Test Suite - Sequential Runner"
  log_msg INFO "======================================"
  log_msg INFO "Configuration:"
  log_msg INFO "  PDB URL: $PDB_URL"
  log_msg INFO "  Reference PDB URL: ${REFERENCE_PDB_URL:-<none>}"
  log_msg INFO "  Number of iterations per mode/type: $NUM_TESTS"
  log_msg INFO "  Base path: $TEST_BASE_PATH"
  log_msg INFO "  Modes: ${TEST_MODES[*]}"
  log_msg INFO "  Types: ${TEST_TYPES[*]}"
  log_msg INFO "  Loop mode: $LOOP_MODE"
  log_msg INFO "======================================"

  local total_test_count=$((NUM_TESTS * ${#TEST_MODES[@]} * ${#TEST_TYPES[@]}))
  if [[ "$LOOP_MODE" == true ]]; then
    log_msg INFO "Sequential test configurations per loop cycle: $total_test_count"
  else
    log_msg INFO "Total test configurations to run: $total_test_count"
  fi

  log_msg INFO "Example folders to be created:"
  for ((i=1; i<=2 && i<=NUM_TESTS; i++)); do
    for mode in "${TEST_MODES[@]}"; do
      for test_type in "${TEST_TYPES[@]}"; do
        if [[ "$LOOP_MODE" == true ]]; then
          log_msg INFO "  - test_${i}_${mode}_${test_type}_YYYYmmdd_HHMMSS/"
        else
          log_msg INFO "  - test_${i}_${mode}_${test_type}/"
        fi
      done
    done
  done
  echo ""
}

run_loop_mode() {
  local cycle=1

  log_msg INFO "Loop mode enabled: tests will run sequentially forever until stopped"
  echo ""

  while true; do
    log_msg INFO "Starting loop cycle #${cycle}"

    for ((i=1; i<=NUM_TESTS; i++)); do
      for mode in "${TEST_MODES[@]}"; do
        for test_type in "${TEST_TYPES[@]}"; do
          log_msg INFO "Loop cycle #${cycle}: running test iteration ${i} for ${mode}/${test_type}"
          run_single_test_sequential "$i" "$mode" "$test_type" || true
          echo ""
        done
      done
    done

    log_msg INFO "Completed loop cycle #${cycle}"
    ((cycle++))
    echo ""
  done
}

main() {
  trap cleanup INT TERM

  parse_arguments "$@"

  print_test_configuration
  check_prerequisites

  mkdir -p "$TEST_BASE_PATH"
  : > "${TEST_BASE_PATH}/global.log"
  log_msg SUCCESS "Test directory created: $TEST_BASE_PATH"

  if [[ "$LOOP_MODE" == true ]]; then
    run_loop_mode
    exit 0
  fi

  log_msg INFO "Running all tests sequentially..."
  echo ""

  local failed_count=0
  for ((i=1; i<=NUM_TESTS; i++)); do
    for mode in "${TEST_MODES[@]}"; do
      for test_type in "${TEST_TYPES[@]}"; do
        run_single_test_sequential "$i" "$mode" "$test_type" || failed_count=$((failed_count + 1))
        echo ""
      done
    done
  done

  log_msg INFO "All tests completed."
  echo ""

  collect_results || failed_count=$?

  local TEST_END_TIME
  TEST_END_TIME=$(date +%s)
  local TOTAL_DURATION=$((TEST_END_TIME - TEST_START_TIME))

  log_msg INFO "======================================"
  log_msg INFO "Total execution time: $(format_duration "$TOTAL_DURATION")"
  log_msg INFO "Results saved to: $TEST_BASE_PATH"
  log_msg INFO "Global log: ${TEST_BASE_PATH}/global.log"
  log_msg INFO "======================================"

  exit $failed_count
}

main "$@"