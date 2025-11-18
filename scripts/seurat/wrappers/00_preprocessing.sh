#!/bin/bash

set -euo pipefail

SMP="$1"

echo "Mapping features of $SMP..." &&
Rscript "$WORK/scripts/seurat/00_map_features.R" 	"$SMP" "$OUT_map/"		"$DIR/$SMP/outs/filtered_feature_bc_matrix" "$TSV" &&
echo "Mapping of $SMP completed✅" &&

echo "Running QC of $SMP..." &&
Rscript "$WORK/scripts/seurat/01_qc.R"			"$SMP" "$OUT_qc/$SMP"		"$OUT_map/${SMP}_mapped.rds" &&
echo "QC of $SMP completed✅" &&

echo "Filtering of $SMP..." &&
Rscript "$WORK/scripts/seurat/02_filter.R"      	"$SMP" "$OUT_filter/$SMP" 	"$OUT_map/${SMP}_mapped.rds" &&
echo "Filtering of $SMP completed✅" &&

echo "Normalizing of $SMP..." &&
Rscript "$WORK/scripts/seurat/03_normalize.R"   	"$SMP" "$OUT_normalize"		"$OUT_filter/$SMP/${SMP}_filtered.rds" &&
echo "Normalization of $SMP completed✅"
