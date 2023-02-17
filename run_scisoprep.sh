#!/usr/bin/env bash
mkdir -p logs

bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -W 24:00 \
  -oo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out \
  -eo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).err \
snakemake \
  -s snake/scisoprep.snake \
  --cores 500 \
  --use-conda \
  --configfile config/config_retina.yaml \
  --profile ~/.config/snakemake/lsf/ \
  -pr \
  --latency-wait 30 \
  --rerun-incomplete \
  "$@"
