#!/bin/sh

for dot in dot_graphs/*.dot; do
  png=${dot%.dot}.png

  echo "Exporting dot graph as png under \"${png}\"..."
  dot -Tpng ${dot} -o ${png}
  echo "Done."
done
