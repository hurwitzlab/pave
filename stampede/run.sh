#!/bin/bash

set -u

function lc() {
   wc -l "$1" | cut -d ' ' -f 1
}

function HELP() {
   printf "Usage:\n  %s -q QUERY_DIR -o OUT_DIR\n\n" $(basename $0)

   echo "Required arguments:"
   echo " -q QUERY_DIR"
   echo " -o OUT_DIR"
   echo ""
   exit 0
}

[[ $# -eq 0 ]] && HELP

QUERY_DIR=""
OUT_DIR="$PWD/pave-out"

while getopts :o:q:h OPT; do
    case $OPT in
      h)
          HELP
          ;;
      o)
          OUT_DIR="$OPTARG"
          ;;
      q)
          QUERY_DIR="$OPTARG"
          ;;
      :)
          echo "Error: Option -$OPTARG requires an argument."
          exit 1
          ;;
      \?)
          echo "Error: Invalid option: -${OPTARG:-""}"
          exit 1
    esac
done

if [[ ! -d $QUERY_DIR ]]; then
    echo QUERY_DIR \"$QUERY_DIR\" does not exist.
    exit 1
fi

[[ ! -d $OUT_DIR ]] && mkdir -p "$OUT_DIR"

DIST_DIR="$OUT_DIR/dist"
QUERY_SKETCH_DIR="$OUT_DIR/sketches"
REPORT_DIR="$OUT_DIR/reports"
REF_SKETCH_DIR="$SCRATCH/imicrobe/mash/sketches"
MASH="$WORK/bin/mash"

if [[ ! -d $REF_SKETCH_DIR ]]; then
  echo REF_SKETCH_DIR \"$REF_SKETCH_DIR\" does not exist.
  exit 1
fi

if [[ ! -d $DIST_DIR ]]; then
  mkdir -p "$DIST_DIR"
fi

if [[ ! -d $QUERY_SKETCH_DIR ]]; then
  mkdir -p "$QUERY_SKETCH_DIR"
fi

if [[ ! -d $REPORT_DIR ]]; then
  mkdir -p "$REPORT_DIR"
fi

#
# Sketch the input files, if necessary
#
ALL_QUERY="$OUT_DIR/all-$(basename $QUERY_DIR)"
if [[ ! -s ${ALL_QUERY}.msh ]]; then
  FILES=$(mktemp)
  find "$QUERY_DIR" -type f > "$FILES"
  NUM_FILES=$(lc "$FILES")

  if [[ $NUM_FILES -lt 1 ]]; then
    echo No files found in QUERY_DIR \"$QUERY_DIR\"
    exit 1
  fi

  echo Sketching NUM_FILES \"$NUM_FILES\"
  while read FILE; do
    SKETCH_FILE="$QUERY_SKETCH_DIR/$(basename $FILE)"
    if [[ -e "${SKETCH_FILE}.msh" ]]; then
      echo SKETCH_FILE \"$SKETCH_FILE.msh\" exists already.
    else
      $MASH sketch -o "$SKETCH_FILE" "$FILE"
    fi
  done < $FILES

  echo Making ALL_QUERY \"$ALL_QUERY\" 

  QUERY_SKETCHES=$(mktemp)
  find "$QUERY_SKETCH_DIR" -name \*.msh > "$QUERY_SKETCHES"
  $MASH paste -l "$ALL_QUERY" "$QUERY_SKETCHES"

  rm "$FILES"
  rm "$QUERY_SKETCHES"
fi
ALL_QUERY=${ALL_QUERY}.msh

#
# The reference genomes ought to have been sketched already
#
ALL_REF="$(dirname $REF_SKETCH_DIR)/all-imicrobe"

if [[ ! -s "${ALL_REF}.msh" ]]; then
  MSH_FILES=$(mktemp)
  find "$REF_SKETCH_DIR" -type f -name \*.msh > $MSH_FILES
  NUM_MASH=$(lc "$MSH_FILES")

  if [[ $NUM_MASH -lt 1 ]]; then
    echo "Found no files in \"$REF_SKETCH_DIR\""
    exit 1
  fi

  echo "Pasting \"$NUM_MASH\" files to ALL_REF \"$ALL_REF\""
  $MASH paste -l "$ALL_REF" "$MSH_FILES"
  rm "$MSH_FILES"
fi
ALL_REF=${ALL_REF}.msh

echo "DIST $(basename $ALL_QUERY) $(basename $ALL_REF)"
DISTANCE_MATRIX="${DIST_DIR}/dist.txt"
$MASH dist -t "$ALL_QUERY" "$ALL_REF" > $DISTANCE_MATRIX

echo "Fixing dist output \"$DIST_DIR\""
./process-dist.pl6 --dist-file=$DISTANCE_MATRIX --out-dir=$OUT_DIR

echo "Done."
