# Genetic Algorism for toxicity prediction

This repository optimizes a set of compounds using a genetic algorithm.  
Additionally, we curate a dataset of compounds used in the validation studies of pharmaceutical safety tests by [JaCVAM](https://www.jacvam.jp/) and compare them with the compounds derived from the genetic algorithm.

## Publication
Not yet.  
  
## Organization
------------
      .
      ├── README.md
      ├── data
      │   ├── processed
      │   │   ├── 0701                      →
      │   │   ├── other
      │   │   ├── tox21
      │   │   └── toxcast
      │   ├── raw
      │   │   ├── Acute_Oral_Toxicity.xlsx
      │   │   ├── DART.xlsx
      │   │   └── validation_test_cas
      │   │       └── 0701.csv
      │   └── result
      │       └── 0701
      │           └── GA_241219_1
      │               ├── GA_chems.png
      │               ├── data.tsv
      │               ├── memo.txt
      │               ├── pareto.png
      │               ├── pareto_front.tsv
      │               ├── phisical_property.png
      │               ├── physical_property_comp_random.png
      │               ├── physical_property_comp_random_only_val.png
      │               ├── score.tsv
      │               ├── shotest_vs_mae.png
      │               ├── shotest_vs_mse.png
      │               ├── shotest_vs_r2.png
      │               ├── structure.png
      │               ├── structure_comp_random.png
      │               ├── structure_comp_random_only_val.png
      │               ├── top_5.tsv
      │               ├── toxicity.png
      │               ├── toxicity_comp_random.png
      │               ├── toxicity_comp_random_only_val.png
      │               ├── toxpred_42
      │               │   ├── best_params.txt
      │               │   ├── best_params_ga.txt
      │               │   ├── difference.png
      │               │   ├── difference_ga.png
      │               │   ├── eval.tsv
      │               │   ├── ga.tsv
      │               │   ├── result.csv
      │               │   ├── result_ga.csv
      │               │   ├── test.tsv
      │               │   └── train.tsv
      │               └── valid_chems.png
      ├── notebook
      │   ├── experiment
      │   │   ├── 0_GA.ipynb
      │   │   ├── 0_dataset_unbalance.ipynb
      │   │   ├── 1_UMAP.ipynb
      │   │   ├── 2_random_generation.ipynb
      │   │   ├── 3_distance_from_pareto.ipynb
      │   │   └── 4_tox_prediction.ipynb
      │   └── preprocess
      │       ├── 0_preprocess_ice.ipynb
      │       ├── 1_preprocess_toxcast.ipynb
      │       ├── 2_preprocess_validation_test.ipynb
      │       ├── 3_preprocess_for_lookup.ipynb
      │       ├── 4_preprocess_for_tox_prediction.ipynb
      │       └── README.md
      ├── requirements.txt
      └── src
         ├── __init__.py
         ├── ga.py
         ├── prep.py
         └── util.py
------------

## Authors
- [Yohei Ohto](https://github.com/YoheiOhto)  
   - main contributor  
- [Tadahaya Mizuno](https://github.com/tadahayamiz)  
  - correspondence  

## Contact
If you have any questions or comments, please feel free to create an issue on github here, or email us:
- oy826c60[at]gmail.com  
- tadahaya[at]gmail.com  
- tadahaya[at]mol.f.u-tokyo.ac.jp  