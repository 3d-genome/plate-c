# 4. A–B difference

**Definition.** log2(AA / BB), the relative strength of A–A versus B–B interactions.

AA and BB come from the same per-replicate saddle analysis as attribute 3 (`../03_compartment_strength`). The `relative_strength_extent10_log2` column of the `AABBABBAvalues_*.tsv` table is log2(AA/BB).

```bash
python common/compare_to_vehicle.py --input AABBABBAvalues_<exp>.tsv \
    --treatment-col drug --vehicle DMSO \
    --value-cols relative_strength_extent10_log2 --output ab_difference_vs_vehicle.tsv
```
