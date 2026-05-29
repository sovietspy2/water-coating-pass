#!/bin/bash
set -euo pipefail

################################################################################
# PASS Test Suite - Parallel Test Runner
# 
# This script runs pass.sh multiple times in parallel with different modes and
# configurations to stress-test the pipeline.
#
# Usage: ./pass_test.sh [OPTIONS]
#
# Options:
#   -u, --url URL           PDB URL to download (default: http://files.rcsb.org/download/1PSV.pdb)
#   -n, --num-tests NUM     Number of parallel tests to run (default: 10)
#   -p, --path PATH         Base path for test directories (default: ./test_runs)
#   -m, --mode MODE         Mode to test: tinker, gromacs, or all (default: all)
#   -t, --type TYPE         Type to test: SHORT, LONG, or all (default: all)
#   -h, --help              Show this help message
#
# Example:
#   ./pass_test.sh -u http://files.rcsb.org/download/1PSV.pdb -n 5 -p /tmp/tests -m all -t all
#
################################################################################

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
readonly PASS_SCRIPT="$PROJECT_ROOT/pipeline/pass.sh"

# Default values
PDB_URL="http://files.rcsb.org/download/1PSV.pdb"
NUM_TESTS=10
TEST_BASE_PATH="./test_runs"
TEST_MODES=("tinker" "gromacs")
TEST_TYPES=("SHORT" "LONG")

# Runtime tracking
declare -a background_pids
declare -A pid_to_test
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
  -u, --url URL           PDB URL to download
                          (default: http://files.rcsb.org/download/1PSV.pdb)
  -n, --num-tests NUM     Number of parallel tests to run
                          (default: 10)
  -p, --path PATH         Base path for test directories
                          (default: ./test_runs)
  -m, --mode MODE         Mode to test: tinker, gromacs, or all
                          (default: all)
  -t, --type TYPE         Type to test: SHORT, LONG, or all
                          (default: all)
  -h, --help              Show this help message

Examples:
  # Default run with 10 tests
  ./pass_test.sh

  # Run 5 tests with gromacs mode only
  ./pass_test.sh -n 5 -m gromacs

  # Run SHORT tests only in /tmp/test_dir
  ./pass_test.sh -p /tmp/test_dir -t SHORT

Notes:
  With -n 2 -t SHORT:
    - 2 tests × 2 modes (tinker, gromacs) = 4 total test runs
    - Each will create a separate folder with its own PDB download
EOF
  exit 0
}

log_msg() {
  local level="$1"
  shift
  local msg="$*"
  local timestamp=$(date '+%Y-%m-%d %H:%M:%S')

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

check_prerequisites() {
  log_msg INFO "Checking prerequisites..."

  if [[ ! -f "$PASS_SCRIPT" ]]; then
    log_msg ERROR "pass.sh not found at: $PASS_SCRIPT"
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
          TEST_TYPES=("SHORT" "LONG")
        else
          TEST_TYPES=("$2")
        fi
        shift 2
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

  if [[ ! -f "$output" ]]; then
    return 1
  fi

  return 0
}

extract_pdb_id() {
  local url="$1"
  basename "$url" .pdb
}

create_test_directory() {
  local test_num="$1"
  local mode="$2"
  local test_type="$3"
  local test_dir="${TEST_BASE_PATH}/test_${test_num}_${mode}_${test_type}"

  mkdir -p "$test_dir"
  echo "$test_dir"
}

run_single_test() {
  local test_num="$1"
  local mode="$2"
  local test_type="$3"
  local pdb_id=$(extract_pdb_id "$PDB_URL")

  local test_dir=$(create_test_directory "$test_num" "$mode" "$test_type")
  local pdb_file="${test_dir}/${pdb_id}.pdb"
  local test_log="${test_dir}/test.log"
  local result_file="${test_dir}/result.txt"

  {
    {
      log_msg INFO "Starting test #${test_num}: MODE=$mode, TYPE=$test_type in $test_dir"

      # Download PDB file
      if ! download_pdb "$PDB_URL" "$pdb_file"; then
        log_msg ERROR "Test #${test_num}: Failed to download PDB"
        echo "FAILED: Download error" > "$result_file"
        exit 1
      fi

      log_msg INFO "Test #${test_num}: PDB downloaded, starting pass.sh"

      # Run pass.sh
      local test_start=$(date +%s)
      if "$PASS_SCRIPT" "$pdb_file" "$mode" "$test_type" >> "$test_log" 2>&1; then
        local test_end=$(date +%s)
        local test_duration=$((test_end - test_start))
        log_msg SUCCESS "Test #${test_num} ($mode/$test_type): Completed in ${test_duration}s"
        echo "SUCCESS: Completed in ${test_duration}s" > "$result_file"
        exit 0
      else
        local exit_code=$?
        log_msg ERROR "Test #${test_num} ($mode/$test_type): Failed with exit code $exit_code"
        echo "FAILED: Exit code $exit_code" > "$result_file"
        tail -20 "$test_log" | log_msg WARN "Last 20 lines of test log:"
        exit 1
      fi
    } 2>&1 | tee -a "${TEST_BASE_PATH}/global.log"
  } &

  # Store the PID and test info
  local pid=$!
  background_pids+=("$pid")
  pid_to_test["$pid"]="test_${test_num}_${mode}_${test_type}"

  echo "$pid"
}

wait_for_tests() {
  local total_tests=${#background_pids[@]}
  local completed=0
  local failed=0

  log_msg INFO "Waiting for $total_tests tests to complete..."

  for pid in "${background_pids[@]}"; do
    if wait "$pid" 2>/dev/null; then
      ((completed++))
    else
      ((failed++))
      ((completed++))
    fi
    local percent=$((completed * 100 / total_tests))
    printf "\r${BLUE}Progress: %d/%d tests completed (%d%%)${NC}" "$completed" "$total_tests" "$percent"
  done

  echo -e "\n"
  return $failed
}

collect_results() {
  log_msg INFO "Collecting test results..."
  echo ""

  local total=0
  local passed=0
  local failed=0
  declare -a failed_tests

  for test_dir in "${TEST_BASE_PATH}"/test_*/; do
    if [[ ! -d "$test_dir" ]]; then
      continue
    fi

    ((total++))
    local test_name=$(basename "$test_dir")
    local result_file="${test_dir}/result.txt"

    if [[ -f "$result_file" ]]; then
      local result=$(cat "$result_file")
      if [[ "$result" =~ "SUCCESS" ]]; then
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

  # Print summary
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
  log_msg INFO "PASS Test Suite - Parallel Runner"
  log_msg INFO "======================================"
  log_msg INFO "Configuration:"
  log_msg INFO "  PDB URL: $PDB_URL"
  log_msg INFO "  Number of iterations per mode/type: $NUM_TESTS"
  log_msg INFO "  Base path: $TEST_BASE_PATH"
  log_msg INFO "  Modes: ${TEST_MODES[*]}"
  log_msg INFO "  Types: ${TEST_TYPES[*]}"
  log_msg INFO "======================================"

  # Calculate total tests
  local total_test_count=$((NUM_TESTS * ${#TEST_MODES[@]} * ${#TEST_TYPES[@]}))
  log_msg INFO "Total test configurations to run: $total_test_count"
  log_msg INFO "Example folders to be created:"
  for ((i=1; i<=2 && i<=NUM_TESTS; i++)); do
    for mode in "${TEST_MODES[@]}"; do
      for test_type in "${TEST_TYPES[@]}"; do
        log_msg INFO "  - test_${i}_${mode}_${test_type}/"
      done
    done
  done
  echo ""
}

main() {
  parse_arguments "$@"

  print_test_configuration

  check_prerequisites

  # Create base test directory and global log
  mkdir -p "$TEST_BASE_PATH"
  : > "${TEST_BASE_PATH}/global.log"
  log_msg SUCCESS "Test directory created: $TEST_BASE_PATH"

  # Launch all tests in parallel
  log_msg INFO "Launching all tests in parallel..."
  echo ""

  for ((i=1; i<=NUM_TESTS; i++)); do
    for mode in "${TEST_MODES[@]}"; do
      for test_type in "${TEST_TYPES[@]}"; do
        run_single_test "$i" "$mode" "$test_type"
        # Small delay to prevent too many simultaneous downloads
        sleep 0.1
      done
    done
  done

  log_msg INFO "All tests launched. Waiting for completion..."
  echo ""

  # Wait for all tests to complete
  local failed_count=0
  if ! wait_for_tests; then
    failed_count=$?
  fi

  # Collect and display results
  if ! collect_results; then
    failed_count=$?
  fi

  local TEST_END_TIME=$(date +%s)
  local TOTAL_DURATION=$((TEST_END_TIME - TEST_START_TIME))

  log_msg INFO "======================================"
  log_msg INFO "Total execution time: ${TOTAL_DURATION}s"
  log_msg INFO "Results saved to: $TEST_BASE_PATH"
  log_msg INFO "Global log: ${TEST_BASE_PATH}/global.log"
  log_msg INFO "======================================"

  exit $failed_count
}

# Run main function
main "$@"