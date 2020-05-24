#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H

#include <vecmath.h>

struct SPHParticle {
    Vector3f position;
    Vector3f velocity;
    Vector3f velocityAtHalfTimeStep;
    Vector3f acceleration;
    double speedOfSound;
    double mass;
    double density;
    double pressure;
    bool isHalfTimeStepVelocityInitialized;
    Vector3f color;
    bool isObstacle;
    bool isVisible = true;
};

#endif
