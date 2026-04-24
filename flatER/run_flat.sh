#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -L)"
MACRO_NAME="Flat_TREEs.C"
FILELIST_DIR="${SCRIPT_DIR}/filelists"

SYSTEM="${1:?missing SYSTEM}"
KIND="${2:?missing KIND}"
TREE="${3:?missing TREE}"
PARTICLE="${4:-}"
PVSNP="${5:-}"

FILELIST_SUBDIR="filelists/${KIND}"
CASETAG="${PARTICLE}${PVSNP}"
FINAL_OUTPUT="flat_${TREE}_${SYSTEM}_${KIND}${CASETAG}.root"
TMP_OUTPUT="$(mktemp "/tmp/flat_${TREE}_${SYSTEM}_${KIND}${CASETAG}.XXXXXX.root")"
INDEXED_OUTPUTS=()
MATCHED_LISTS=()
SYSTEM_MATCH="${SYSTEM,,}"
KIND_MATCH="${KIND,,}"
TREE_MATCH="${TREE,,}"
PARTICLE_MATCH="${PARTICLE,,}"
PVSNP_MATCH="${PVSNP,,}"

# For DATA, only the system matters. For MC ntmix, particle matters too. For other MC, tree is enough.
while IFS= read -r list_path; do
  list_name="$(basename "${list_path}")"
  list_match="${list_name,,}"

  if [[ "${KIND}" == "DATA" ]]; then
    if [[ "${list_match}" == *"${KIND_MATCH}"* && "${list_match}" == *"${SYSTEM_MATCH}"* ]]; then
      MATCHED_LISTS+=("${list_path}")
    fi
  else
    if [[ "${TREE}" == "ntmix" ]]; then
      if [[ "${list_match}" == *"${KIND_MATCH}"* && "${list_match}" == *"${SYSTEM_MATCH}"* && "${list_match}" == *"${TREE_MATCH}"* && "${list_match}" == *"${PARTICLE_MATCH}"* ]]; then
        if [[ "${PVSNP_MATCH}" == *"nonprompt"* ]]; then
          [[ "${list_match}" == *"nonprompt"* ]] && MATCHED_LISTS+=("${list_path}")
        else
          [[ "${list_match}" != *"nonprompt"* ]] && MATCHED_LISTS+=("${list_path}")
        fi
      fi
    else
      if [[ "${list_match}" == *"${KIND_MATCH}"* && "${list_match}" == *"${SYSTEM_MATCH}"* && "${list_match}" == *"${TREE_MATCH}"* ]]; then
        MATCHED_LISTS+=("${list_path}")
      fi
    fi
  fi
done < <(find "${FILELIST_SUBDIR}" -type f -name '*.txt' | sort)

cd "${SCRIPT_DIR}"

# Run the flattener once per matched list and tag each chunk with _0, _1, _2, ...
for idx in "${!MATCHED_LISTS[@]}"; do
  FILELIST="${MATCHED_LISTS[$idx]}"
  NUN="_${idx}"
  OUTPUT_CHUNK="flat_${TREE}_${SYSTEM}_${KIND}${CASETAG}${NUN}.root"

  root -l -b -q "${MACRO_NAME}(\"${FILELIST}\",\"${NUN}\",\"${TREE}\",\"${SYSTEM}\",\"${KIND}\",\"${PARTICLE}\",\"${PVSNP}\")"
  INDEXED_OUTPUTS+=("${OUTPUT_CHUNK}")
done

# Merge all chunk outputs into one final file, then remove the indexed temporary files.
rm -f "${FINAL_OUTPUT}"
if [[ "${#INDEXED_OUTPUTS[@]}" -eq 1 ]]; then
  mv -f "${INDEXED_OUTPUTS[0]}" "${FINAL_OUTPUT}"
else
  hadd -f "${TMP_OUTPUT}" "${INDEXED_OUTPUTS[@]}"
  mv -f "${TMP_OUTPUT}" "${FINAL_OUTPUT}"
  rm -f "${INDEXED_OUTPUTS[@]}"
fi
