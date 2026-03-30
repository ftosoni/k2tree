#!/bin/bash
# test_consistency.sh
# Automates the cross-checking of different k2tree formats

set -e

BUILD_DIR=${1:-"."}
ROOT_DIR=${2:-".."}

# Define paths to executables
# On Windows/MinGW, suffix might be .x.exe
K2SPARSE="$BUILD_DIR/k2sparse.x"
K2CPDF="$BUILD_DIR/k2cpdf.x"
B128SPARSE="$BUILD_DIR/b128sparse.x"
MATRIXCMP="$BUILD_DIR/matrixcmp.x"

# Compatibility check for Windows (.exe suffix)
if [[ ! -f "$K2SPARSE" && -f "$K2SPARSE.exe" ]]; then
    K2SPARSE="$K2SPARSE.exe"
    K2CPDF="$K2CPDF.exe"
    B128SPARSE="$B128SPARSE.exe"
    MATRIXCMP="$MATRIXCMP.exe"
fi

# Create a temporary directory for test files
TEMP_DIR="tests/temp_consistency"
mkdir -p "$TEMP_DIR"

INPUT_TXT="$TEMP_DIR/input.txt"
TEST_K2="$TEMP_DIR/test.k2"
TEST_CK2="$TEMP_DIR/test.ck2"
TEST_B128="$TEMP_DIR/test.b128"
DECODED_TXT="$TEMP_DIR/decoded.txt"

echo "--- Generating test matrix (16x16 sparse) ---"
cat <<EOF > "$INPUT_TXT"
0 0
1 1
2 2
3 3
4 4
5 5
6 6
7 7
8 8
9 9
10 10
11 11
12 12
13 13
14 14
15 15
0 15
15 0
5 10
EOF

echo -e "\n--- Testing Standard K2-TREE (.k2) ---"
"$K2SPARSE" "$INPUT_TXT" -o "$TEST_K2"
"$K2SPARSE" -d "$TEST_K2" -o "$DECODED_TXT"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: K2-TREE matches input."

echo -e "\n--- Testing Compressed K2-DFS (.ck2) ---"
"$K2CPDF" -o "$TEST_CK2" "$TEST_K2"
BACKP="${TEST_CK2}.p"
"$K2SPARSE" -d "$TEST_CK2" -I "$BACKP" -o "$DECODED_TXT"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: K2-DFS (Compressed) matches input."

echo -e "\n--- Testing B128 Baseline (.b128) ---"
"$B128SPARSE" "$INPUT_TXT" -o "$TEST_B128"
"$B128SPARSE" -d "$TEST_B128" -o "$DECODED_TXT"
"$MATRIXCMP" "$INPUT_TXT" "$DECODED_TXT"
echo "SUCCESS: B128 matches input."

echo -e "\nALL FORMATS ARE CONSISTENT!"

# Cleanup
# rm -rf "$TEMP_DIR"
