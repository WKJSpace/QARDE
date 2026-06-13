open_component /tmp/qarde_workspace_csim_top_mm -reset -flow_target vivado
set_top qarde_accel_mm
add_files /home/hlab/wei/workstation/workspace/QARDE/src/qarde_accel_mm.cpp
add_files -tb /home/hlab/wei/workstation/workspace/QARDE/testbench/tb_top_mm.cpp
set_part xczu48dr-ffvg1517-2-e
create_clock -period 2.8
csim_design
exit
