#!/bin/bash

# main suite versions

create_versions () {
    # the first parameter is a list, passed by name reference
    local -n list=$1 # the caller passes in the array to fill
    var_base="ec=x0.5_mtz=0_u=0.8"
    var_1=("ml=x0.2" "ml=x0.5" "ml=x0.7")
    var_2=("lc=x0" "lc=x0.5" "lc=x1")
    var_3=("uc=x0" "uc=x0.2" "uc=x0.8")
    var_4=("rec=70" "rec=160" "rec=250")  
    for a in ${var_1[@]}; do 
        for b in ${var_2[@]}; do 
            for c in ${var_3[@]}; do 
                for d in ${var_4[@]}; do
                    list+=("${a} ${b} ${c} ${d} ${var_base}");
                done
            done
        done
    done
}