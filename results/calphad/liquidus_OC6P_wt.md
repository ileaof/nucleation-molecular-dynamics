# Al-Nb — OC6P (OpenCalphad 6.083), base ALFENB-B2-4SL.TDB, composições em wt.%

Macro: `calphad_oc/alnb_wt.OCM` (rodar no WSL Ubuntu: `~/OC6HOME/oc6P alnb_wt.OCM`).
Método: a fase que aparece é fixada com quantidade 0 (`set status phase <fase>=fix 0`) e T é liberada (`set cond t=none`).

| Nb na liga | T_L estável (Al₃Nb aparece) | Peritético L + Al₃Nb → FCC_A1 | Nb no líquido no peritético | Nb na FCC_A1 nascente | T_L metaestável FCC_A1 (só LIQUID + FCC) | Nb na FCC nesse ponto |
|---|---|---|---|---|---|---|
| 5,00 wt.% | 1539,27 K | **934,51 K** | **0,0062 wt.%** | 0,516 wt.% | 1059,66 K | 21,3 wt.% |
| 5,16 wt.% | 1544,20 K | 934,51 K | 0,0062 wt.% | 0,516 wt.% | 1061,69 K | 21,5 wt.% |
| 11,01 wt.% | 1675,20 K | 934,51 K | 0,0062 wt.% | 0,516 wt.% | 1122,16 K | 27,0 wt.% |

Al₃Nb: 53,4 wt.% Nb em todos os casos. DTA (script de nucleação, L218): T_L = 936,58 K.
Os valores coincidem com o pycalphad (`scripts/calphad_liquidus.py`, base 2SL).
