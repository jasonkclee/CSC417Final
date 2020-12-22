# Fluid simulation
Uses LibIgl for rendering. Make sure to have the libigl installed in the root directory, as folder "libigl", or clone from the repo with `git clone --recursive`
Link to video: https://drive.google.com/file/d/1RpY4kBs0ZHhLSZXkEahUdkOaYhZUQQk1/view?usp=sharing
Link to repo: https://github.com/jasonkclee/CSC417Final
### Compilation
`mk dir build`
`cd build`
`cmake .. -DCMAKE_BUILD_TYPE=Release`
`make`
Run with `./fluid`
Use the `t` and `r` key to spawn more particles
### Optional Arguments
`-s [grid size]`
`-n [number of particles]`
`-s [grid size]`
`-dt [delta time]`
`-g [gravity]`
`-pic` Uses the PIC method instead of FLIP
EX: `./fluid -s 10 -n 3`