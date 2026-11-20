# Capstone: Which nuclear features distinguish malignant from benign breast masses?

## Biological question
Fine-needle-aspirate images of breast masses describe the cell nuclei (size, shape, texture).
**Which of these 30 nuclear features differ between malignant and benign masses, and by how much?**

## Data
Wisconsin Diagnostic Breast Cancer (WDBC) dataset, UCI Machine Learning Repository (public).
569 samples (212 malignant, 357 benign), 30 numeric features (mean, standard error and "worst" of 10 measurements).
Loaded through `scikit-learn`, so no download step is needed.

## Methods
1. Each feature compared between groups with a two-sided **Mann-Whitney U test** (no normality assumption).
2. Effect size: **rank-biserial correlation** (positive = higher in malignant, range -1 to 1).
3. **Multiple-testing control** over 30 tests: Bonferroni (family-wise error) and **Benjamini-Hochberg FDR** (q < 0.05).
4. Feature-feature Spearman correlations computed to check redundancy.

## How to run
```
pip install -r requirements.txt
python analysis.py
```
Outputs go to `results/`: `feature_tests.csv`, `feature_correlations.csv`, `volcano.png`, `top_feature_boxplot.png`.
(Or open in Google Colab, upload `analysis.py`, and run `!pip install -r requirements.txt` then `!python analysis.py`.)

## Results
- Raw p < 0.05: 27/30 features. After Bonferroni: 26/30. After BH-FDR (q < 0.05): 27/30.
- Strongest differences: *worst perimeter, worst radius, worst area, worst concave points, mean concave points*
  (rank-biserial about 0.93-0.95; all higher in malignant masses).
- Not significant after correction: *mean fractal dimension, texture error, smoothness error*.
- Multiple-testing correction barely changed the result because the effects are very large; it mainly removed one borderline feature (Bonferroni).

## Interpretation (cautious)
Malignant masses show larger nuclei with more irregular outlines (more concave points), consistent with known
cytological features of malignancy such as nuclear enlargement and pleomorphism. Size features separate the groups strongly.
Shape and size are not independent findings, though, so this should not be read as 27 separate biological signals.

## Limitations and confounders
- **Correlated features:** radius, perimeter and area are geometrically linked (median |r| between features = 0.42, and much higher within size features). The 30 tests are not independent, so the number of "significant" features overstates independent signals. BH assumes independence or positive dependence; it is still reasonable here, but counts remain inflated.
- **Association, not causation:** a statistical difference does not explain why nuclei enlarge.
- **No clinical covariates:** age, menopausal status, mass type, imaging device and operator are not in the data and could confound differences.
- **Single source/cohort:** one institution's samples, so results may not generalise. No external validation was done.
- **Sample selection:** these are aspirates already referred for diagnosis, so the benign/malignant balance (37% malignant) does not reflect screening populations.
- **Statistical vs practical significance:** with extreme effects, p-values (down to ~1e-80) are not informative beyond "very small". Effect sizes are more meaningful.
- Possible batch effects (measurement day/slide preparation) cannot be assessed.
