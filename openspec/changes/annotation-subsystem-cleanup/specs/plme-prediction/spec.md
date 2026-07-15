## REMOVED Requirements

### Requirement: PLM-Effector off by default
**Reason**: PLM-Effector is removed from ssign entirely. It over-called ~25% of a proteome and has been off at every tier since 2026-06-17; carrying dead-but-loadable code works against the low-maintenance release goal.
**Migration**: No action needed for default runs (PLM-Effector was already off by default, so default output is unchanged). The `--skip-plm-effector` / `--no-skip-plm-effector` and `--plme-batch-size` flags are removed; passing them now produces an unknown-argument error.

### Requirement: PLM-Effector remains installable
**Reason**: PLM-Effector is removed entirely; its package (`plm_effector/`), wrapper (`run_plm_effector.py`), merge helper, weights references, and dependency declarations are deleted.
**Migration**: The `secretion-classifier-dataset` project must source any PLM-Effector-derived feature independently of ssign (full removal chosen deliberately on 2026-07-15).

### Requirement: PLM-Effector positivity confidence gate
**Reason**: PLM-Effector is removed, so there is no PLM-Effector positive/negative call to gate.
**Migration**: None. Secretion voting in `cross_validate` now counts DeepLocPro + DeepSecE only (`n_prediction_tools_agreeing` ranges 0-2); SignalP remains evidence-only; `is_secreted` is still "at least one predictor flags".
