open_component /tmp/qarde_workspace_csim_ems -reset -flow_target vivado
set_top main
add_files -tb /home/hlab/wei/workstation/workspace/QARDE/testbench/tb_ems_algorithm.cpp
set_part xczu48dr-ffvg1517-2-e
create_clock -period 2.8
csim_design
exit
