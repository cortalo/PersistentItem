## Persistent Item


### About this repo

This repo contains the algorithms of finding persistent items using cold filter framework.

| Task                      | Algorithms                               |
| ------------------------- | ---------------------------------------- |
| Finding Persistent Item | PIE |


### Requirements

- The gather-and-report part of CF use SIMD instructions to achieve high speed, so the cpu must support SSE2 instruction set.
- cmake >= 2.6
- g++ (MSVC is not supported currently.)

### How to build

The project is tested on Microsoft Visual Studio 2017. 

Change the location of dataset in main.cpp, 424 line. Then you can run the program.

### Related paper
[Cold Filter: A Meta-Framework for Faster and More Accurate Stream Processing](http://www.zhouy.me/paper/cf-sigmod18.pdf)
