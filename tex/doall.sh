#! /bin/sh

for file in ToF_beam LandauProfiles ; do
  pdflatex ${file}.tex

done
