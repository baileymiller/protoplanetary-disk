# protoplanetary-disk
A tool for rendering protoplanetary disks. Uses PBRT to create the renderings and python to generate the voxel grids.


## Building

```
	mkdir build
	cd build
	cmake ..
```

## Building
From within the python directory, call the following:

```
python create_disk.py
```

Within `create_disk.py` you can adjust the parameters of the disk or render a range of parameter values.


![alt text](./edge_on_disk.png)