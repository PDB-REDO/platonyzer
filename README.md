platonyzer
==========

This is the repository for platonyzer, an application to create restraints for metal sites. This program is part of the [PDB-REDO](https://pdb.redo.eu/) suite
of programs.

Installation
------------

To install, first install [`libcifpp`](https://github.com/PDB-REDO/libcifpp) and [`libpdb-redo`](https://github.com/PDB-REDO/libpdb-redo) then use

```bash
git clone https://github.com/PDB-REDO/platonyzer.git
cd platonyzer
mkdir build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
sudo cmake --install build
```

