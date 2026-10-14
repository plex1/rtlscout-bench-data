#!/bin/bash
# usage: ./sim.sh design.sv
D=${1:-design.sv}
B=${D%.sv}
which verilator >/dev/null 2>&1 && {
  rm -rf obj_$B
  verilator --binary --timing -Wno-fatal -j 4 --top tb -o tb_$B tb.sv $D -Mdir obj_$B > ver_$B.log 2>&1
  ./obj_$B/tb_$B 2>&1 | tail -20
} || echo "no verilator"
