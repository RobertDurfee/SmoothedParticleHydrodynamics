#ifndef SPH_OBSTACLE_H
#define SPH_OBSTACLE_H

#include <vector>
#include "vecmath.h"
#include "sphparticle.h"

struct SPHObstacle {
    Vector3f position;
    std::vector<SPHParticle*> particles;
    bool isVisible = true;
};

#endif
