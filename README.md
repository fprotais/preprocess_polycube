# pre-process with naive flagging for polycube

Pre-process a mesh so that the interior may not have impossible configuration for a flagging. I give a naive flagging with it, so it is possible to obtain a hexmesh with my other code: [https://github.com/fprotais/polycube_withHexEx](https://github.com/fprotais/polycube_withHexEx).

Modification done: no tet should have more than 1 boundary facet and if an edge is not on boundary, then both of its point can't be on the boundary. Split are made accordingly.

The code is very slow because I am lazy. I might improve it later. 

# Use CMake to build the project:
```sh
git clone --recurse-submodules https://github.com/fprotais/preprocess_polycube &&
cd preprocess_polycube &&
mkdir build &&
cd build &&
cmake -DCMAKE_BUILD_TYPE=Release .. &&
make -j 
```

# Running the code :

```sh
./preprocess ../S1.mesh  res.mesh  res.flags(optional)
```
For the supported mesh formats, see [ultimaille](https://github.com/ssloy/ultimaille). 
