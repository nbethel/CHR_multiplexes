#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH -t 0:10:00
#SBATCH --gres=gpu:titan:1
#SBATCH -c 3
#SBATCH --output=mpnn_run.out

source activate mlfold
export PATH=/home/nbethel/mpnn-master/:$PATH
python /home/nbethel/mpnn-master/mpnn_run_tied.py \
        --max_length 10000 \
        --checkpoint_path '/projects/ml/struc2seq/data_for_complexes/training_scripts/paper_experiments/model_outputs/p10/checkpoints/epoch51_step255000.pt' \
        --hidden_dim 192 \
        --num_layers 3 \
        --protein_features 'full' \
        --jsonl_path='pdbs.jsonl' \
        --chain_id_jsonl 'pdbs_masked.jsonl'  \
        --fixed_positions_jsonl '' \
        --out_folder='output' \
        --num_seq_per_target 4 \
        --sampling_temp="0.1" \
        --batch_size 4 \
        --omit_AAs 'X' \
        --backbone_noise 0.05 \
        --decoding_order 'random' \
        --bias_AA_jsonl '' \
        --num_connections 64 \
        --tied_positions_jsonl 'pdbs_tied.jsonl'

