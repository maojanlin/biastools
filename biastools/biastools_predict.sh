path_out=$1
sample_id=$2
run_id=$3
flag_real=$4
report_real=$5
report_simulation=$6
prefix=${path_out}/${sample_id}


if [[ ${report_real} == 'none' ]]; then
    report_real=${prefix}.real.${run_id}.bias
fi

if [[ ${flag_real} == 1 || ${report_simulation} == 'none' ]]; then
    echo "[Biastools] Real report bias prediction."
    python3 predict_model.py -rr ${report_real} -out ${prefix}.real.${run_id}
else
    echo "[Biastools] Bias prediction based on simulation report!"
    python3 predict_experiment.py -sr ${report_simulation} \
                                  -rr ${report_real} \
                                  -out ${prefix}.sim.${run_id}
fi