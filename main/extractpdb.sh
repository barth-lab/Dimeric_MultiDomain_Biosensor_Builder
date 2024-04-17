#!/bin/bash

# extract all pdbs from a silent file
while getopts ":R:S:" opt; do
  case $opt in
    R)
      R=$OPTARG # rosetta location
      ;;
    S)
      S=$OPTARG # silent file input
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

$R/source/bin/extract_pdbs.linuxgccrelease -in:file:silent_struct_type binary -in:file:silent $S
