#!/usr/bin/env bash
bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -W 24:00 \
  -oo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).out \
  -eo logs/snakelog.$(date +%Y-%m-%d.%H-%M-%S).err \
snakemake \
  -s snake/scisoprep.snake \
  --configfile config/config.yaml \
  --profile ~/.config/snakemake/lsf/ \
  --restart-times 3 \
  -pr \
  --latency-wait 30 \
  --rerun-incomplete \
  "$@"
