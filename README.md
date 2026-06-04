The repository is provided as supplementary material for the MA thesis 
"Longitudinal Trajectories of Turn-Taking in Japanese Child-Caregiver Dyads".

## RL_all_annotations.csv
The dataset contains processed and anonymized response-latency observations used in the analyses. 
No raw audio recordings are included.
For RQ1, observations should be filtered to the response-latency range [-4, 7]. 
For RQ2, observations with missing values in either `PrevSelfLatency` or `PrevInterLatency` 
should be excluded before model fitting.

## 2.2.2.2 Computation of response latency
This script computes response latency (RL_s) from turn-level annotations. For each recording, turns are ordered by onset time, and response latency is calculated as the difference between the onset of the current turn and the offset of the immediately preceding turn by the other speaker. Positive values indicate gaps, while negative values indicate overlap.
