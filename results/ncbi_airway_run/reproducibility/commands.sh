python3 de_pipeline.py \
  --counts /home/yer_kanat/Downloads/lol/depipeline1/data/ncbi_airway_counts.csv \
  --metadata /home/yer_kanat/Downloads/lol/depipeline1/data/ncbi_airway_metadata.csv \
  --formula "~ batch + condition" \
  --contrast "condition,treated,control" \
  --output /home/yer_kanat/Downloads/lol/depipeline1/results/ncbi_airway_run
