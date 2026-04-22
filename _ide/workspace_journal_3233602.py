# 2026-03-26T22:53:24.831647276
import vitis

client = vitis.create_client()
client.set_workspace(path="QARDE")

cfg = client.get_config_file(path="/home/asap_jupiter/wei/workspace/QARDE/QARDE/hls_config.cfg")

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp", "../src/qarde_tools.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_config.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_config.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_ems_cnu.hpp", "../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_tools.hpp", "../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_gf.hpp", "../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/nb_ldpc.hpp", "../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde.hpp", "../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_dec.hpp", "../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_accel.cpp"])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_value(section="hls", key="syn.top", value="fp_ems_accel")

comp = client.get_component(name="QARDE")
comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

cfg.set_value(section="hls", key="syn.top", value="")

cfg.set_value(section="hls", key="syn.top", value="qarde_accel")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

comp.run(operation="SYNTHESIS")

vitis.dispose()

