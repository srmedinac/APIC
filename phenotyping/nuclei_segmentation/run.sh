#!/bin/sh
for i in /wsi/subfolders/with/tiles/*; do
  # Check if it's a directory
  if [ -d "$i" ]; then
    # Extract just the folder name without the full path
    folder_name=$(basename "$i")

    # Define output directory
    out_dir="output/dir/path/$folder_name"

    # Check if output directory already exists, if so, skip
    if [ -d "$out_dir" ]; then
      echo "Output directory $out_dir already exists, skipping..."
      continue # Skip the rest of this iteration
    fi

    # If output directory does not exist, create it
    mkdir -p "$out_dir"

    python src/nucleusSegmentationTiles.py "/path/to/histoprep/tiles/$folder_name/tiles" '.jpeg' 'models/hovernet_fast_pannuke_type_tf2pytorch.tar' 'fast' 6 '/path/to/types/type_info_pannuke.json' "$out_dir"
  fi
done