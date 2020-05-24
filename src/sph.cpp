#include "sph.h"

#include "vertexrecorder.h"
#include "util.h"

SPH::SPH() {
    // Create a bunch of particles with random positions (all zero velocity).
    std::vector<Vector3f> positions;
    for (int i = 0; i < this->numberOfParticles; i++) {
        float xPosition = rand_uniform(this->xMinOpeningPosition, this->xMaxOpeningPosition);
        float yPosition = rand_uniform(this->yMinOpeningPosition, this->yMaxOpeningPosition);
        float zPosition = rand_uniform(this->zMinOpeningPosition, this->zMaxOpeningPosition);
        positions.push_back(Vector3f(xPosition, yPosition, zPosition));
    }
    this->addFluidParticles(positions);
    // Create a dam to demonstrate fluid qualities.
    this->addDamObstacle();
    // Add boundary particles if requested (didn't seem to work super well).
    if (this->enableBoundaryParticles) {
        this->addBoundaryParticles();
    }
}

void SPH::addDamObstacle() {
    // Calculate number of particles in wall (y-z plane).
    int zParticles = (int)ceil((this->zMaxBoundary - this->zMinBoundary) / (2 * this->particleRadius));
    int yParticles = (int)ceil((this->yMaxBoundary - this->zMinBoundary) / (2 * this->particleRadius));
    // Should have some thickness otherwise particles push through.
    int wallWidth = 3;
    std::vector<Vector3f> positions;
    for (int xIndex = 0; xIndex < wallWidth; xIndex++) {
        // Calculate x position (particles should be flush side-by-side).
        double xPosition = 2 * xIndex * this->particleRadius;
        for (int zIndex = 0; zIndex <= zParticles; zIndex++) {
            // Calculate z position (particles should be flush side-by-side).
            double zPosition = 2 * zIndex * this->particleRadius + this->zMinBoundary;
            for (int yIndex = 0; yIndex <= yParticles; yIndex++) {
                double yPosition = 2 * yIndex * this->particleRadius + this->yMinBoundary;
                // Don't add if within the dam opening.
                if (zPosition >= this->zMinDamOpening && zPosition <= this->zMaxDamOpening && yPosition >= this->yMinDamOpening && yPosition <= this->yMaxDamOpening) {
                    continue;
                }
                positions.push_back(Vector3f(xPosition, yPosition, zPosition));
            }
        }
    }
    this->addObstacleParticles(positions);
}

void SPH::addBoundaryParticles() {
    int xParticles = (int)ceil((this->yMaxBoundary - this->zMinBoundary) / (2 * this->particleRadius));
    int yParticles = (int)ceil((this->yMaxBoundary - this->zMinBoundary) / (2 * this->particleRadius));
    int zParticles = (int)ceil((this->zMaxBoundary - this->zMinBoundary) / (2 * this->particleRadius));
    // Left Wall (y-z plane)
    std::vector<Vector3f> positions;
    for (int yIndex = 0; yIndex <= yParticles; yIndex++) {
        // Calculate y position (particles should be flush side-by-side).
        double yPosition = 2 * yIndex * this->particleRadius + this->yMinBoundary;
        for (int zIndex = 0; zIndex <= zParticles; zIndex++) {
            // Calculate z position (particles should be flush side-by-side).
            double zPosition = 2 * zIndex * this->particleRadius + this->zMinBoundary;
            positions.push_back(Vector3f(this->xMinBoundary, yPosition, zPosition));
        }
    }
    int obstacleID = this->obstacles.size();
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
    obstacleID++;
    // Right Wall (y-z plane).
    positions.clear();
    for (int yIndex = 0; yIndex <= yParticles; yIndex++) {
        // Calculate y position (particles should be flush side-by-side).
        double yPosition = 2 * yIndex * this->particleRadius + this->yMinBoundary;
        for (int zIndex = 0; zIndex <= zParticles; zIndex++) {
            // Calculate z position (particles should be flush side-by-side).
            double zPosition = 2 * zIndex * this->particleRadius + this->zMinBoundary;
            positions.push_back(Vector3f(this->xMaxBoundary, yPosition, zPosition));
        }
    }
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
    obstacleID++;
    // Floor (x-z plane)
    positions.clear();
    for (int xIndex = 0; xIndex <= xParticles; xIndex++) {
        // Calculate x position (particles should be flush side-by-side).
        double xPosition = 2 * xIndex * this->particleRadius + this->xMinBoundary;
        for (int zIndex = 0; zIndex <= zParticles; zIndex++) {
            // Calculate z position (particles should be flush side-by-side).
            double zPosition = 2 * zIndex * this->particleRadius + this->zMinBoundary;
            positions.push_back(Vector3f(xPosition, this->yMinBoundary, zPosition));
        }
    }
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
    obstacleID++;
    // Ceiling (x-z plane).
    positions.clear();
    for (int xIndex = 0; xIndex <= xParticles; xIndex++) {
        // Calculate x position (particles should be flush side-by-side).
        double xPosition = 2 * xIndex * this->particleRadius + this->xMinBoundary;
        for (int zIndex = 0; zIndex <= zParticles; zIndex++) {
            // Calculate z position (particles should be flush side-by-side).
            double zPosition = 2 * zIndex * this->particleRadius + this->zMinBoundary;
            positions.push_back(Vector3f(xPosition, this->yMaxBoundary, zPosition));
        }
    }
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
    obstacleID++;
    // Back Wall (x-y plane).
    positions.clear();
    for (int xIndex = 0; xIndex <= xParticles; xIndex++) {
        // Calculate x position (particles should be flush side-by-side).
        double xPosition = 2 * xIndex * this->particleRadius + this->xMinBoundary;
        for (int yIndex = 0; yIndex <= yParticles; yIndex++) {
            // Calculate y position (particles should be flush side-by-side).
            double yPosition = 2 * yIndex * this->particleRadius + this->yMinBoundary;
            positions.push_back(Vector3f(xPosition, yPosition, this->zMinBoundary));
        }
    }
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
    obstacleID++;
    // Front Wall (x-y plane).
    positions.clear();
    for (int xIndex = 0; xIndex <= xParticles; xIndex++) {
        // Calculate x position (particles should be flush side-by-side).
        double xPosition = 2 * xIndex * this->particleRadius + this->xMinBoundary;
        for (int yIndex = 0; yIndex <= yParticles; yIndex++) {
            // Calculate y position (particles should be flush side-by-side).
            double yPosition = 2 * yIndex * this->particleRadius + this->yMinBoundary;
            positions.push_back(Vector3f(xPosition, yPosition, this->zMaxBoundary));
        }
    }
    this->addObstacleParticles(positions);
    // Make sure to hide these particles from view.
    for (unsigned int i = 0; i < this->obstacles[obstacleID]->particles.size(); i++) {
        this->obstacles[obstacleID]->particles[i]->isVisible = false;
    }
}

void SPH::update(double dt) {
    double timeLeft = dt;
    while (timeLeft > 0.0) {
        // Get current density, pressure, and acceleration (must be done in this order).
        this->updateFluidDensityAndPressure();
        this->updateFluidAcceleration();
        // We don't necessarily want to take the full time-step `dt`.
        double timeStep = this->calculateTimeStep();
        timeLeft -= timeStep;
        if (timeLeft < 0.0) {
            timeStep = timeStep + timeLeft;
            timeLeft = 0.0;
        }
        // Use the acceleration to update position (step by the adjusted time step).
        this->updateFluidPosition(timeStep);
    }
}

void SPH::draw(GLProgram& glProgram) {
    SPHParticle *particle;
    // We want to draw all particles (including obstacles).
    for (unsigned int i = 0; i < this->allParticles.size(); i++) {
        particle = this->allParticles[i];
        // Skip particles marked as invisible (this would be the boundary particles--if used).
        if (!particle->isVisible) {
            continue;
        }
        // Particles might have difference colors (e.g. obstacle vs. fluid).
        glProgram.updateMaterial(particle->color);
        glProgram.updateModelMatrix(Matrix4f::translation(particle->position));
        // Draw a sphere representing the particle.
        drawSphere(this->particleRadius, this->particleSlices, this->particleStacks);
    }
}

void SPH::addFluidParticles(std::vector<Vector3f> positions) {
    for (unsigned int i = 0; i < positions.size(); i++) {
        // All particles should start with no velocity.
        this->addFluidParticle(positions[i], Vector3f::ZERO);
    }
}

void SPH::addObstacleParticles(std::vector<Vector3f> positions) {
    SPHObstacle *obstacle = new SPHObstacle();
    SPHParticle *particle;
    for (unsigned int i = 0; i < positions.size(); i++) {
        // Add obstacle particle to simulation and save to current obstacle.
        particle = this->addObstacleParticle(positions[i]);
        obstacle->particles.push_back(particle);
    }
    // Save obstacle (this isn't absolutely necessary, but might be helpful in the future).
    this->obstacles.push_back(obstacle);
}

void SPH::addFluidParticle(Vector3f position, Vector3f velocity) {
    SPHParticle *particle = this->createSPHParticle(position, velocity);
    // Obstacles are treated as stationary, be sure to not treat these particles as obstacles.
    particle->isObstacle = false;
    // Add to both fluid and all particles (helpful to differentiate when updating).
    this->fluidParticles.push_back(particle);
    this->allParticles.push_back(particle);
}

SPHParticle* SPH::addObstacleParticle(Vector3f position) {
    SPHParticle *particle = createSPHObstacleParticle(position);
    // Make sure particle is obstacle so it remains fixed.
    particle->isObstacle = true;
    // Special color to highlight.
    particle->color = this->obstacleColor;
    // Add to both obstacle and all particles (helpful to differentiate when updating).
    this->obstacleParticles.push_back(particle);
    this->allParticles.push_back(particle);
    return particle;
}

SPHParticle* SPH::createSPHParticle(Vector3f position, Vector3f velocity) {
    SPHParticle *particle = new SPHParticle();
    // Initialize with shared defaults among all particles (e.g. mass, density, etc.).
    particle->position = position;
    particle->velocity = velocity;
    particle->velocityAtHalfTimeStep = Vector3f::ZERO;
    particle->isHalfTimeStepVelocityInitialized = false;
    particle->speedOfSound = 0;
    particle->acceleration = Vector3f::ZERO;
    particle->density = this->initialDensity;
    particle->mass = this->particleMass;
    particle->pressure = 0.;
    particle->color = this->particleColor;
    return particle;
}

SPHParticle* SPH::createSPHObstacleParticle(Vector3f position) {
    // Same as fluid particle, essentially.
    SPHParticle *particle = this->createSPHParticle(position, Vector3f::ZERO);
    return particle;
}

void SPH::updateFluidDensityAndPressure() {
    SPHParticle *iParticle, *jParticle;
    Vector3f radius;
    // Loop over *all* particles (including obstacles because we use neighboring 
    //   obstacle particle pressures/densities when calculating fluid acceleration).
    for (unsigned int i = 0; i < this->allParticles.size(); i++) {
        iParticle = this->allParticles[i];
        double density = 0.0;
        // Again, loop over *all* particles.
        for (unsigned int j = 0; j < this->allParticles.size(); j++) {
            jParticle = this->allParticles[j];
            radius = iParticle->position - jParticle->position;
            double distance = radius.abs();
            // Only skip the particle if distance is 0 (i.e. comparing the same
            //   particle with itself) or the particle is outside the smoothing
            //   radius (if don't ignore, things blow up).
            if (distance == 0.0 || distance > this->smoothingRadius) {
                continue;
            }
            double distanceSquared = Vector3f::dot(radius, radius);
            double difference = this->smoothingRadiusSquared - distanceSquared;
            density += jParticle->mass * this->poly6Coefficient * pow(difference, 3.);
        }
        // If the density is lower than the initial, we will get negative 
        //   pressures (this will cause issues down the line) so take max.
        iParticle->density = fmax(density, this->initialDensity);
        iParticle->pressure = this->pressureCoefficient * (iParticle->density - this->initialDensity);
    }
}

void SPH::updateFluidAcceleration() {
    // Make sure to call `updateFluidDensityAndPressure()` before calling this.
    SPHParticle *iParticle, *jParticle;
    Vector3f acceleration;
    Vector3f radius;
    Vector3f velocityDifference;
    // Only consider the fluid particles, obstacle particles have no acceleration.
    for (unsigned int i = 0; i < this->fluidParticles.size(); i++) {
        iParticle = this->fluidParticles[i];
        acceleration = Vector3f::ZERO;
        // But, consider neighbors from *all* particles as they will have an impact
        //   on the current particle's motion.
        for (unsigned int j = 0; j < this->allParticles.size(); j++) {
            jParticle = this->allParticles[j];
            radius = iParticle->position - jParticle->position;
            double distance = radius.abs();
            // Only skip the particle if distance is 0 (i.e. comparing the same
            //   particle with itself) or the particle is outside the smoothing
            //   radius (if don't ignore, things blow up).
            if (distance == 0.0 || distance > this->smoothingRadius) {
                continue;
            }
            float inverseDistance = 1 / distance;
            radius = inverseDistance * radius;
            float difference = this->smoothingRadius - distance;
            float spikey = this->spikeyGradientCoefficient * pow(difference, 2.);
            float massRatio = jParticle->mass / iParticle->mass;
            float p = (iParticle->pressure + jParticle->pressure) / (2 * iParticle->density * jParticle->density);
            acceleration -= this->accelerationCoefficient * massRatio * p * spikey * radius;
            // There is only a viscous force between two fluid particles, not 
            //   between a fluid particle and an obstacle.
            if (!jParticle->isObstacle) {
                float laplacian = this->viscosityLaplacianCoefficient * difference;
                velocityDifference = jParticle->velocity - iParticle->velocity;
                acceleration += this->viscosityCoefficient * massRatio * (1 / jParticle->density) * laplacian * velocityDifference;
            }
        }
        acceleration += this->gravity;
        // We allow the boundary to have some repulsive force on the particles to
        //   try to prevent particles from sticking to walls.
        acceleration += this->calculateBoundaryAcceleration(iParticle);
        double accelerationMagnitude = acceleration.abs();
        // Apply damping to slow particle motion if configured to do so
        Vector3f damping = this->motionDampingCoefficient * iParticle->velocity;
        if (damping.abs() > accelerationMagnitude) {
            acceleration = Vector3f::ZERO;
        } else {
            acceleration -= damping;
        }
        // If acceleration is too large, things will quickly become unstable,
        //   this should prevent unstable simulations
        if (accelerationMagnitude > this->maximumAcceleration) {
            acceleration = (acceleration / accelerationMagnitude) * this->maximumAcceleration;
        }
        iParticle->acceleration = acceleration;
    }
}

void SPH::updateFluidPosition(double dt) {
    // Make sure to call `updateFluidAcceleration()` before calling this.
    SPHParticle *particle;
    // Only consider fluid particles (as obstacles don't move)
    for (unsigned int i = 0; i < this->fluidParticles.size(); i++) {
        particle = this->fluidParticles[i];
        // We do a form of trapezoidal integration, make sure halfway derivative
        //   is initialized.
        if (particle->isHalfTimeStepVelocityInitialized) {
            particle->velocityAtHalfTimeStep += dt * particle->acceleration;
        } else {
            particle->velocityAtHalfTimeStep = particle->velocity + (0.5 * dt) * particle->acceleration;
            particle->isHalfTimeStepVelocityInitialized = true;
        }
        // Update the position and velocity for time step
        particle->position += dt * particle->velocityAtHalfTimeStep;
        particle->velocity += particle->velocityAtHalfTimeStep + (0.5 * dt) * particle->acceleration;
        // Make sure the velocity doesn't get too large otherwise things might
        //   become unstable. This attempts to prevent that.
        if (particle->velocity.abs() > this->maximumVelocity) {
            Vector3f unitVelocity = particle->velocity.normalized();
            particle->velocity = this->maximumVelocity * unitVelocity;
        }
        // Make sure particles don't leave the box. (This might be the reason
        //   some particles stick in the edges/corners of the box.)
        this->enforceFluidParticlePositionBounds(particle);
    }
}

double SPH::calculateTimeStep() {
    // Taking too large of a time step can result in unstable simulation. This
    //   method is provided by Ryan L. Guy, "The goal of choosing a time step is 
    //   to choose the largest value such that a particle will move less than one 
    //   smoothing radius in a single time step." (http://rlguy.com/sphfluidsim/)
    double maxVelocitySquared = 0.0;
    double maxSpeedOfSoundSquared = 0.0;
    double maxAccelerationSquared = 0.0;
    SPHParticle *particle;
    for (unsigned int i = 0; i < this->fluidParticles.size(); i++) {
        particle = this->fluidParticles[i];
        double velocitySquared = Vector3f::dot(particle->velocity, particle->velocity);
        double speedOfSoundSquared = this->evaluateSpeedOfSoundSquared(particle);
        double accelerationSquared = Vector3f::dot(particle->acceleration, particle->acceleration);
        if (velocitySquared > maxVelocitySquared) {
            maxVelocitySquared = velocitySquared;
        }
        if (speedOfSoundSquared > maxSpeedOfSoundSquared) {
            maxSpeedOfSoundSquared = speedOfSoundSquared;
        }
        if (accelerationSquared > maxAccelerationSquared) {
            maxAccelerationSquared = accelerationSquared;
        }
    }
    double maxVelocity = sqrt(maxVelocitySquared);
    double maxSpeedOfSound = sqrt(maxSpeedOfSoundSquared);
    double maxAcceleration = sqrt(maxAccelerationSquared);
    double velocityStep = this->courantSafetyFactor * this->smoothingRadius / fmax(1.0, maxVelocity);
    double speedOfSoundStep = this->courantSafetyFactor * this->smoothingRadius / maxSpeedOfSound;
    double accelerationStep = sqrt(this->smoothingRadius / maxAcceleration);
    return fmax(this->minTimeStep, fmin(fmin(velocityStep, speedOfSoundStep), accelerationStep));
}

inline double lerp(double x1, double x2, double t) {
    // Simple linear interpolation
    return x1 + t * (x2 - x1);
}

Vector3f SPH::calculateBoundaryAcceleration(SPHParticle *particle) {
    // This provides a increasing force from the boundary the closer the particle
    //   gets to the boundary (not currently used)
    double radius = this->boundaryForceRadius;
    Vector3f position = particle->position;
    Vector3f acceleration = Vector3f::ZERO;
    if (position.x() < this->xMinBoundary + radius) {
        double distance = fmax(0.0, position.x() - this->xMinBoundary);
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance / radius);
        acceleration += Vector3f(force / particle->mass, 0.0, 0.0);
    } else if (position.x() > this->xMaxBoundary - radius) {
        double distance = fmax(0.0, this->xMaxBoundary - position.x());
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance / radius);
        acceleration += Vector3f(-force / particle->mass, 0.0, 0.0);
    }
    if (position.y() < this->yMinBoundary + radius) {
        double distance = fmax(0.0, position.y() - this->yMinBoundary);
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance / radius);
        acceleration += Vector3f(0.0, force / particle->mass, 0.0);
    } else if (position.y() > this->yMaxBoundary - radius) {
        double distance = fmax(0.0, this->yMaxBoundary - position.y());
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance/radius);
        acceleration += Vector3f(0.0, -force / particle->mass, 0.0);
    }
    if (position.z() < this->zMinBoundary + radius) {
        double distance = fmax(0.0, position.z() - this->zMinBoundary);
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance / radius);
        acceleration += Vector3f(0.0, 0.0, force / particle->mass);
    } else if (position.z() > this->zMaxBoundary - radius) {
        double distance = fmax(0.0, this->zMaxBoundary - position.z());
        double force = lerp(this->maxBoundaryForce, this->minBoundaryForce, distance / radius);
        acceleration += Vector3f(0.0, 0.0, -force / particle->mass);
    }
    return acceleration;
}

void SPH::enforceFluidParticlePositionBounds(SPHParticle *particle) {
    // Just make sure the particles don't get outside of the bounding box.
    double epsilon = 0.001;
    float d = this->boundaryDampingCoefficient;
    if (particle->position.x() < this->xMinBoundary) {
        particle->position = Vector3f(this->xMinBoundary + epsilon, particle->position.y(), particle->position.z());
        particle->velocity = Vector3f(-d * particle->velocity.x(), particle->velocity.y(), particle->velocity.z());
    } else if (particle->position.x() > this->xMaxBoundary) {
        particle->position = Vector3f(this->xMaxBoundary - epsilon, particle->position.y(), particle->position.z());
        particle->velocity = Vector3f(-d * particle->velocity.x(), particle->velocity.y(), particle->velocity.z());
    }
    if (particle->position.y() < this->yMinBoundary) {
        particle->position = Vector3f(particle->position.x(), this->yMinBoundary + epsilon, particle->position.z());
        particle->velocity = Vector3f(particle->velocity.x(), -d * particle->velocity.y(), particle->velocity.z());
    } else if (particle->position.y() > this->yMaxBoundary) {
        particle->position = Vector3f(particle->position.x(), this->yMaxBoundary - epsilon, particle->position.z());
        particle->velocity = Vector3f(particle->velocity.x(), -d * particle->velocity.y(), particle->velocity.z());
    }
    if (particle->position.z() < this->zMinBoundary) {
        particle->position = Vector3f(particle->position.x(), particle->position.y(), this->zMinBoundary + epsilon);
        particle->velocity = Vector3f(particle->velocity.x(), particle->velocity.y(), -d * particle->velocity.z());
    } else if (particle->position.z() > this->zMaxBoundary) {
        particle->position = Vector3f(particle->position.x(), particle->position.y(), this->zMaxBoundary - epsilon);
        particle->velocity = Vector3f(particle->velocity.x(), particle->velocity.y(), -d * particle->velocity.z());
    }
}

double SPH::evaluateSpeedOfSoundSquared(SPHParticle *particle) {
    // Make sure this calculation doesn't blow up.
    if (particle->density < 0.00001) {
        return 0.0;
    }
    return this->ratioOfSpecificHeats * particle->pressure / particle->density;
}
