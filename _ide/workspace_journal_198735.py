# 2026-03-26T22:21:38.146692260
import vitis

client = vitis.create_client()
client.set_workspace(path="QARDE")

comp = client.create_hls_component(name = "QARDE",cfg_file = ["hls_config.cfg"],template = "empty_hls_component")

cfg = client.get_config_file(path="/home/asap_jupiter/wei/workspace/QARDE/QARDE/hls_config.cfg")

cfg.set_values(key="syn.file", values=["../src/qarde_init.hpp"])

cfg.set_values(key="syn.file", values=[])

cfg.set_values(key="syn.blackbox.file", values=[])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp"])

cfg.set_values(key="syn.file", values=["../src/qarde_vnu.hpp", "../src/qarde_init.hpp", "../src/qarde_perm.hpp", "../src/qarde_ems_cnu.hpp", "../src/quarde_dec.hpp"])

vitis.dispose()

