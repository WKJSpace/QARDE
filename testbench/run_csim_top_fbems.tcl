open_component /tmp/qarde_workspace_csim_top_fbems -reset -flow_target vivado
set_top qarde_accel_fbems
add_files /home/hlab/wei/workstation/workspace/QARDE/src/qarde_accel_fbems.cpp
add_files -tb /home/hlab/wei/workstation/workspace/QARDE/testbench/tb_top_fbems.cpp
set_part xczu48dr-ffvg1517-2-e
create_clock -period 2.8
csim_design
exit
