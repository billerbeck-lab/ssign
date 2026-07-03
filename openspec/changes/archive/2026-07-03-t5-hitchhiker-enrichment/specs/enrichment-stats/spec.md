## MODIFIED Requirements

### Requirement: Per-type aggregation with autotransporter self-detection
The enrichment test SHALL pool a genome's systems of each SS type. For window types (T1SS, T2SS, T3SS, T4SS, pT4SSt, T5bSS, T6SSi, T6SSii) the observed signal SHALL be secreted-predicted proteins within the configured gene window (default +/-3) of any component of that type. For autotransporter types (T5aSS, T5cSS) the test SHALL emit **two** results distinguished by the `mode` column: a `self` result (the component itself is predicted outer-membrane-or-extracellular for DLP, or secreted-typed for DSE, because the autotransporter is both machinery and substrate) AND a `window` "hitchhiker" result (secreted-predicted proteins within the configured gene window of that type's components, testing the hypothesis that neighbours piggyback through the T5 pore). The hitchhiker window SHALL use the same window mask, predictors, and combined rule as other window types. DSE SHALL NOT be tested for T3SS.

#### Scenario: Window type uses neighborhood
- **WHEN** the tested SS type is a window type
- **THEN** the observed count SHALL be secreted-predicted proteins within the configured gene window of that type's components

#### Scenario: Autotransporter type emits self and hitchhiker results
- **WHEN** the tested SS type is T5aSS or T5cSS
- **THEN** the enrichment output SHALL contain a `mode=self` result (component self-detection) AND a `mode=window` hitchhiker result (secreted-predicted proteins within the configured gene window of that type's components), for each applicable predictor

#### Scenario: Hitchhiker combined track uses the window rule
- **WHEN** the combined (union) track is computed for a T5aSS/T5cSS hitchhiker (`mode=window`) result
- **THEN** it SHALL be the DLP-or-DSE union, the same rule as other window types (whereas the `mode=self` combined track remains the DLP-or-SignalP union)

#### Scenario: T3SS excluded from DSE
- **WHEN** the predictor is DSE
- **THEN** no enrichment row SHALL be produced for T3SS
