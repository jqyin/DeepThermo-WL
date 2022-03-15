#ifndef MAIN_H
#define MAIN_H

#ifdef CRUSHER
#define GPUperNode 8
#else
#define GPUperNode 6
#endif

enum SamplingMode {metropolis, WLprior, WLdos, WLproduction};

#endif
