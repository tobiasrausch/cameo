#!/bin/bash

# Download Jaspar motifs
curl --output jaspar.txt 'https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt'
gzip jaspar.txt
