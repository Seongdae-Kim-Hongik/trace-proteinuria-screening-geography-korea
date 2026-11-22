#!/bin/bash
# run_std3_chain.sh - Full model chain for the principal version (std3) and the sensitivity versions (s577, s2, s8, sc, sp): panels (80), optional checks,
#  81 (principal model and multiverse), 81R, 85, 86, 87, 88g, 88h, 88k, 89g. 88S and 88V are run separately afterwards.
#  Requires ST_GRAPH (principal adjacency graph). R jobs run through rslot.sh with one INLA thread; return codes are written to run_std3_rc.log.
set -x
S="${ST_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}"; cd "$S"; RS="$S/01_스크립트/rslot.sh"; export SK_RSLOTS=7
T=NBKID; GR="${ST_GRAPH:?set ST_GRAPH to the principal adjacency graph file}"; D=261122
RCLOG="$S/04_분석/$T/run_std3_rc.log"; rcl(){ echo "$(date '+%F %T') rc=$2 $1" >> "$RCLOG"; }
PIDS=(); NAMES=(); bgadd(){ PIDS+=("$1"); NAMES+=("$2"); }
waitall(){ local i rc; for i in "${!PIDS[@]}"; do wait "${PIDS[$i]}"; rc=$?; rcl "${NAMES[$i]}" "$rc"; done; PIDS=(); NAMES=(); }
wl(){ while :; do l=$(sysctl -n vm.loadavg | awk '{print $2}'); awk -v l="$l" 'BEGIN{exit !(l<14)}' && break; echo "[wait] load1 $l"; sleep 60; done; }
q(){ while [ "$(jobs -rp | wc -l)" -ge 3 ]; do sleep 10; done; wl; }
L="$S/04_분석/$T"; LG="${TMPDIR:-/tmp}/nbkid_logs"; mkdir -p "$LG"
ST_VER=std3 python3 01_스크립트/80_build_disease_panel_v.py $T > $LG/run_80.log 2>&1; rc=$?; rcl 80_build_panel $rc; [ $rc -eq 0 ] || { echo "80 FAILED"; exit 1; }
[ -z "${ST_GATE_DIR:-}" ] || python3 "$ST_GATE_DIR/data_gate.py" $L/panel_NBKID_std3.csv --region fbd_region > $LG/data_gate_NBKID_std3_261121.txt 2>&1; rcl data_gate $?
[ -z "${ST_GATE_DIR:-}" ] || python3 "$ST_GATE_DIR/panel_gate.py" $L/panel_NBKID_std3.csv --name-col fbd_region --cols population,E,cases,offset --positive population,E,offset --years 2011-2024 > $LG/panel_gate_NBKID_std3_261121.txt 2>&1; rcl panel_gate $?
python3 - "$L/panel_INFO_NBKID_std3.json" <<'PYEOF' || { rcl vif_gate 1; echo "VIF GATE FAILED"; exit 1; }
import json,sys; J=json.load(open(sys.argv[1])); o=J.get("vif_over_threshold",[]); print("vif_over_threshold",o); sys.exit(1 if o else 0)
PYEOF
rcl vif_gate 0
export ST_VER=std3 ST_THREADS=1 ST_GRAPH="$GR"
wl; (cd 01_스크립트 && $RS Rscript 81_fit_inla.R $T) > $LG/run_std3_81.log 2>&1 & P81=$!
sleep 30; wl; (cd 01_스크립트 && $RS Rscript 81R_nbkid_refit_binom.R $T) > $LG/run_std3_81R.log 2>&1 & bgadd $! 81R
( for V in s577 s2 s8 sc sp; do ST_VER=$V python3 01_스크립트/80_build_disease_panel_v.py $T && wl && (cd 01_스크립트 && ST_VER=$V ST_THREADS=1 ST_MV_N=2 ST_GRAPH=$([ $V = s577 ] || echo "$GR") $RS Rscript 81_fit_inla.R $T); done > $LG/run_std3_sensver.log 2>&1 ) & bgadd $! sensver
wait $P81; rc=$?; rcl 81_full $rc; [ $rc -eq 0 ] || { echo "81 FAILED"; exit 1; }
echo "81 done"
python3 01_스크립트/87_window_grid.py $T std3 > $LG/run_std3_87.log 2>&1; rcl 87_window_grid $?
wl; (cd 01_스크립트 && TAG=$T ST_VER=std3 ST_GRAPH=$GR $RS Rscript 85_window_maps.R $T std3) > $LG/run_std3_85.log 2>&1; rcl 85_window_maps $?
python3 01_스크립트/86_yearly_maps_ehsa.py $T std3 > $LG/run_std3_86.log 2>&1; rcl 86_yearly_maps_ehsa $?
python3 01_스크립트/88g_sens_prep.py $T std3 KIDS,NBBP > $LG/run_std3_88gprep.log 2>&1; rc=$?; rcl 88g_sens_prep $rc; [ $rc -eq 0 ] || exit 1
cd 01_스크립트; export TAG=$T NCO=KIDS,NBBP
for G in 1 2 3 4 5 6; do q; GROUP=$G $RS Rscript 88g_sens_suite.R > $LG/run_std3_88g_$G.log 2>&1 & bgadd $! 88g_$G; sleep 5; done
q; $RS Rscript 88h_ppc_varpart.R > $LG/run_std3_88h.log 2>&1 & bgadd $! 88h
for G in nb struct review conf; do q; $RS Rscript 88k_struct_kh.R $G > $LG/run_std3_88k_$G.log 2>&1 & bgadd $! 88k_$G; sleep 5; done
q; $RS Rscript 88k_struct_kh_anem.R kh > $LG/run_std3_88k_kh.log 2>&1 & bgadd $! 88k_kh
waitall
if grep -q 'KH_typeII FAILED\|KH_typeIV FAILED' $LG/run_std3_88k_kh.log; then KH_INIT=0 $RS Rscript 88k_struct_kh_anem.R kh > $LG/run_std3_88k_kh_init0.log 2>&1; rcl 88k_kh_init0 $?; fi
cd ..; TAG=$T ST_VER=std3 python3 01_스크립트/89g_sens_combine_nbkid.py $T std3 > $LG/run_std3_89g.log 2>&1; rc=$?; rcl 89g_sens_combine_nbkid $rc
cd 01_스크립트; TAG=$T ST_VER=std3 ST_GRAPH=$GR $RS Rscript 88V_nbkid_extra.R $T > $LG/run_std3_88V.log 2>&1; rcl 88V_extra $?
TAG=$T ST_VER=std3 ST_GRAPH=$GR $RS Rscript 88S_nbkid_reverse_posthoc.R > $LG/run_std3_88S.log 2>&1; rcl 88S_posthoc $?; cd ..
for f in "$LG"/*; do cp -X "$f" "$L/"; done
echo "ALL rc=$rc"; grep -l "FAILED" $LG/run_std3_88*.log; true
