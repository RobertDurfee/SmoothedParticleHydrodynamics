#ifndef SPH_H
#define SPH_H

#include <vector>
#include <vecmath.h>
#include "sphparticle.h"
#include "sphobstacle.h"
#include "glprogram.h"

const double PI = 3.14159265358979323846;

class SPH {

    // This implementation is based on the work done by Ryan L. Guy 
    //   (http://rlguy.com/sphfluidsim/).
    
public:

    SPH();

    void update(double dt);
    void draw(GLProgram& glProgram);

    void addFluidParticles(std::vector<Vector3f> positions);
    void addFluidParticle(Vector3f position, Vector3f velocity);

    void addObstacleParticles(std::vector<Vector3f> positions);
    SPHParticle* addObstacleParticle(Vector3f position);

    void addDamObstacle();
    void addBoundaryParticles();

private:

    SPHParticle* createSPHParticle(Vector3f position, Vector3f velocity);
    SPHParticle* createSPHObstacleParticle(Vector3f position);

    void updateFluidDensityAndPressure();
    void updateFluidAcceleration();
    void updateFluidPosition(double dt);

    double calculateTimeStep();
    Vector3f calculateBoundaryAcceleration(SPHParticle *particle);
    void enforceFluidParticlePositionBounds(SPHParticle *particle);
    double evaluateSpeedOfSoundSquared(SPHParticle *particle);

    // Simulation Constants
    const int numberOfParticles = 500;
    const double smoothingRadius = 0.2;
    const double smoothingRadiusSquared = pow(this->smoothingRadius, 2);
    const Vector3f gravity = Vector3f(0.0f, -9.8f, 0.0f);
    const double courantSafetyFactor = 1.0;
    const double minTimeStep = 1.0 / 240.0;
    const double initialDensity = 50.0;
    const double pressureCoefficient = 20.0;
    const double particleMass = 1.0;
    const double motionDampingCoefficient = 0.0;
    const double boundaryDampingCoefficient = 0.6;
    const double ratioOfSpecificHeats = 1.0;
    const double accelerationCoefficient = 0.3;
    const double viscosityCoefficient = 0.0018;
    const double maximumVelocity = 75.0;
    const double maximumAcceleration = 75.0;
    const bool enableBoundaryParticles = false; // This didn't seem to have benefits

    // Kernel Constants (precompute once)
    const double poly6Coefficient = 315.0 / (64.0 * PI * pow(this->smoothingRadius, 9.0));
    const double spikeyGradientCoefficient = -45.0 / (PI * pow(this->smoothingRadius, 6.0));
    const double viscosityLaplacianCoefficient = 45.0 / (PI * pow(this->smoothingRadius, 6.0));

    // Boundary Constraints
    const double boundaryForceRadius = 0.1;
    const double minBoundaryForce = 0.0;
    const double maxBoundaryForce = 0.0;
    const double xMinBoundary = -0.5;
    const double xMaxBoundary = 2.0;
    const double yMinBoundary = -0.5;
    const double yMaxBoundary = 2.0;
    const double zMinBoundary = -0.5;
    const double zMaxBoundary = 0.5;

    // Particle Generator Constraints
    const double xMinOpeningPosition = -0.5;
    const double xMaxOpeningPosition = 0.0;
    const double yMinOpeningPosition = -0.5;
    const double yMaxOpeningPosition = 0.5;
    const double zMinOpeningPosition = -0.5;
    const double zMaxOpeningPosition = 0.5;

    // Dam Obstacle Constraints
    double yMinDamOpening = -0.5;
    double yMaxDamOpening = 0.0;
    double zMinDamOpening = -0.25;
    double zMaxDamOpening = 0.25;

    // Graphics Constants
    const Vector3f particleColor = Vector3f(0.4f, 0.7f, 1.f); // Blueish
    const Vector3f obstacleColor = Vector3f(0.4f, 0.4f, 0.4f); // Grayish
    const double particleRadius = 0.04;
    const int particleSlices = 8;
    const int particleStacks = 8;

    // Obstacles (not really necessary to keep these but makes life a little easier)
    std::vector<SPHObstacle*> obstacles;

    // State
    std::vector<SPHParticle*> fluidParticles; // Moving particles
    std::vector<SPHParticle*> obstacleParticles; // Stationary particles
    std::vector<SPHParticle*> allParticles;
    
};

#endif
