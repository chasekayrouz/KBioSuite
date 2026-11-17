"""Run external engines without importing their dependencies into the desktop host.

Inference sockets are blocked. MSA retrieval is an explicitly separate operation.
This file must remain runnable with the external environment's Python.
"""

import json
import runpy
import socket
import sys
from pathlib import Path


def block_network():
    original_connect = socket.socket.connect
    original_connect_ex = socket.socket.connect_ex

    def connect(instance, address):
        if instance.family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError(
                "Kmatch local inference cannot contact the network. Install assets first."
            )
        return original_connect(instance, address)

    def connect_ex(instance, address):
        if instance.family in (socket.AF_INET, socket.AF_INET6):
            raise RuntimeError("Network access is disabled during Kmatch local inference.")
        return original_connect_ex(instance, address)

    socket.socket.connect = connect
    socket.socket.connect_ex = connect_ex


def main():
    operation = sys.argv.pop(1)
    if operation == "msa":
        sequence_file, output = map(Path, sys.argv[1:3])
        from boltz.data.msa.mmseqs2 import run_mmseqs2

        output.mkdir(parents=True, exist_ok=True)
        # Public service receives the sequence only after explicit GUI/CLI selection.
        result = run_mmseqs2(
            sequence_file.read_text().strip(),
            prefix=str(output / "service"),
            host_url="https://api.colabfold.com",
        )
        a3m = result[0] if isinstance(result, tuple) else result
        partial = output / "query.a3m.part"
        partial.write_text(a3m[0])
        partial.replace(output / "query.a3m")
        (output / "source.json").write_text(
            json.dumps(
                {
                    "mode": "network-assisted-msa",
                    "service": "https://api.colabfold.com",
                    "structure_inference": "local",
                },
                indent=2,
            )
            + "\n"
        )
    else:
        block_network()
        if operation == "boltz":
            from boltz.main import cli

            cli()
        elif operation == "af3":
            source = Path(sys.argv.pop(1))
            sys.argv[0] = str(source / "run_alphafold.py")
            runpy.run_path(sys.argv[0], run_name="__main__")
        else:
            raise ValueError("Unknown engine operation.")


if __name__ == "__main__":
    main()
