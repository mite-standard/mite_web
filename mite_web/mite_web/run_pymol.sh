#!/bin/bash

mkdir -p ./mite_web/static/img

for pdb in ./mite_web/data/pdb/*.pdb; do
  filename=$(basename -- "$pdb");
  filename=${filename%.pdb}
  pymol -c "$pdb" -d "bg_color white; hide everything; show cartoon; spectrum count, red blue; set opaque_background, 0; png ./mite_web/static/img/$filename, 0, 0, -1, ray=1; quit;";
  done
