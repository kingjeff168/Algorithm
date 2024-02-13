#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <time.h>


#define leftWall   -1
#define rightWall  -2
#define upWall     -3
#define downWall   -4


typedef struct item
{
    double particleXPosition;
    double particleYPosition;
    double particleXVelocity;
    double particleYVelocity;
    double timeCollision;
    int countWall;
    int countParticleCollision;
    int index;
}Item;

// collision between two objects and time
typedef struct collision
{
    int object1;
    int object2;
    double timeCollision;
}Collision;

void Particle_Wall_Collision(double radius, Item * p, Collision * c, double time, int length, int width)
{
    double timeX, timeY;

    // Use the formula to calculate the time of collision in x and y direction.
    if (p->particleXVelocity < 0)
    {
        timeX = ((radius - p->particleXPosition) / p->particleXVelocity) + time;
    }
    else
    {
        timeX = ((length - radius - p->particleXPosition) / p->particleXVelocity) + time;
    }

    if (p->particleYVelocity < 0)
    {
        timeY = ((radius - p->particleYPosition) / p->particleYVelocity) + time;
    }
    else
    {
        timeY = ((width - radius - p->particleYPosition) / p->particleYVelocity) + time;
    }

    // Compare timeX, timeY, particleXVelocity, and particleYVelocity to determine the first wall the particle collides.
    if (timeX < timeY)
    {
        if (p->particleXVelocity < 0)
        {
            c->object1 = p->index;
            c->object2 = leftWall;
            c->timeCollision = timeX;
        }
        else
        {
            c->object1 = p->index;
            c->object2 = rightWall;
            c->timeCollision = timeX;
        }
    }
    else
    {
        if (p->particleYVelocity < 0)
        {
            c->object1 = p->index;
            c->object2 = downWall;
            c->timeCollision = timeY;
        }
        else
        {
            c->object1 = p->index;
            c->object2 = upWall;
            c->timeCollision = timeY;
        }
    }
}

void Particles_Collision(Item * p1, Item * p2, Collision * cStruct, double timeReference, double radius)
{
    double timeResult, timeResult1, timeResult2, tempV1, tempV2;
    double a, b, c, d, i, j, k, l;
    a = p1->particleXVelocity - p2->particleXVelocity;
    b = p1->particleXPosition - p2->particleXPosition;
    c = p1->particleYVelocity - p2->particleYVelocity;
    d = p1->particleYPosition - p2->particleYPosition;

    // Calculate the time of collision using the formula answer of quadratic equation.
    i = (a * a) + (c * c);
    j = 2*((a * b) + (c * d));
    k = (b * b) + (d * d) - (4 * radius * radius);
    l = sqrt((j * j) - (4 * i * k));

    if (l >= 0)
    {
        timeResult1 = (((-1) * j) + sqrt((j * j) - (4 * i * k))) / (2 * i);
        timeResult2 = (((-1) * j) - sqrt((j * j) - (4 * i * k))) / (2 * i);
        if ((timeResult1 > 0) && (timeResult1 < timeResult2))
        {
            timeResult = timeResult1 + timeReference;
        }
        else if ((timeResult2 > 0) && (timeResult2 < timeResult1))
        {
            timeResult = timeResult2 + timeReference;
        }
        else
        {
            // Even l >= 0, the time might be negative, so I assign a very large value for timeResult to avoid being sorted to the first one.
            timeResult = 10000;
        }

        cStruct->object1 = p1->index;
        cStruct->object2 = p2->index;
        cStruct->timeCollision = timeResult;
    }
    else
    {
        // l < 0 means no solution, so I assign a very large value for timeResult.
        cStruct->object1 = p1->index;
        cStruct->object2 = p2->index;
        timeResult = 10000;
        cStruct->timeCollision = timeResult;
    }
}

// This function is to update the position, velocity, and time stamp of particle-wall collision.
void particleWallUpdate(Item * p, int index, double time, double radius, int length, int width)
{
    // Update the position and velocity in x direction.
    p->particleXPosition = p->particleXPosition + (p->particleXVelocity * (time - p->timeCollision));

    if ((index == leftWall) || (index == rightWall))
    {
        p->particleXVelocity = (-1) * (p->particleXVelocity);
    }

    // Update the position and velocity in y direction.
    p->particleYPosition = p->particleYPosition + (p->particleYVelocity * (time - p->timeCollision));
            
    if ((index == upWall) || (index == downWall))
    {
        p->particleYVelocity = (-1) * (p->particleYVelocity);
    }

    p->timeCollision = time;

}

// This function is to update the position, velocity, and time stamp of two particles' collision.
void particlesUpdate(Item * p1, Item * p2, double time, double radius, int length, int width)
{
    double tempV1, tempV2;
    // Update the position.
    p1->particleXPosition = p1->particleXPosition + p1->particleXVelocity * (time - p1->timeCollision);
    p1->particleYPosition = p1->particleYPosition + p1->particleYVelocity * (time - p1->timeCollision);
    p2->particleXPosition = p2->particleXPosition + p2->particleXVelocity * (time - p2->timeCollision);
    p2->particleYPosition = p2->particleYPosition + p2->particleYVelocity * (time - p2->timeCollision);

    // Update the velocity.
    tempV1 = p1->particleXVelocity;
    tempV2 = p1->particleYVelocity;
    p1->particleXVelocity = p2->particleXVelocity;
    p1->particleYVelocity = p2->particleYVelocity;
    p2->particleXVelocity = tempV1;
    p2->particleYVelocity = tempV2;

    p1->timeCollision = time;
    p2->timeCollision = time;
}

// This function is to use update the position of the particle at currentTime.
void particleDataUpdate(Item * p, double time)
{
    p->particleXPosition = p->particleXPosition + (p->particleXVelocity * (time - p->timeCollision));
    p->particleYPosition = p->particleYPosition + (p->particleYVelocity * (time - p->timeCollision));
    p->timeCollision = time;
}


void insertionSort(Collision * array, int n) 
{
    Collision key;
    int i, j;
    for (i = 1; i < n; i++) 
    {
        key = array[i];
        j = i - 1;
        while (j >= 0 && array[j].timeCollision > key.timeCollision) 
        {
            array[j + 1] = array[j];
            j = j - 1;
        }
        array[j + 1] = key;
    }
}



int main(int argc, char *argv[])
{
    FILE *file;
    char *filename;
    int num_elements = 0;
    int count, xLength, yWidth, totalCollisionStruct, index, index2;
    double radius, currentTime, endTime;

    filename = argv[1];

    // Open the file.
    file = fopen(filename, "r");

    // Read the number of particles.
    if (fscanf(file, "%d", &count) != 1) 
    {
        printf("Error reading count from the file.\n");
        fclose(file);
        return 1;
    }

    // Read the radius.
    if (fscanf(file, "%lf", &radius) != 1) 
    {
        printf("Error reading radius from the file.\n");
        fclose(file);
        return 1;
    }

    // Read the box's length and width.
    if (fscanf(file, "%d %d", &xLength, &yWidth) != 2) 
    {
        printf("Error reading box's dimension from the file.\n");
        fclose(file);
        return 1;
    }

    // The total amount of struct collision.
    totalCollisionStruct = (count * (count + 1)) / 2;
    Item * particles = (Item *)malloc(count * sizeof(Item));
    Collision * collision = (Collision *)malloc(totalCollisionStruct * sizeof(Collision));


    // Read all items and put into struct array.
    while (fscanf(file, "%lf %lf %lf %lf", &particles[num_elements].particleXPosition, &particles[num_elements].particleYPosition, &particles[num_elements].particleXVelocity, &particles[num_elements].particleYVelocity) != EOF) 
    {
        num_elements++;
    }
    fclose(file);


    // Read end time from command-line arguments.
    endTime = atoi(argv[2]);

    // Initialize the other values in the struct.
    for (int i = 0; i < count; i++)
    {
        particles[i].timeCollision = 0;
        particles[i].countParticleCollision = 0;
        particles[i].countWall = 0;
        particles[i].index = i;
    }

    for (int i = 0; i < totalCollisionStruct; i++)
    {
        collision[i].object1 = 0;
        collision[i].object2 = 0;
        collision[i].timeCollision = 0;
    }
    

    // Conduct the first calculation of the collision.
    // Use "currentTime" to assign the time of the collision after sorting.
    currentTime = 0;
    for (int j = 0; j < count; j++)
    {
        Particle_Wall_Collision(radius, &particles[j], &collision[j], currentTime, xLength, yWidth);
    }

    // Because the start position of struct for Particles_Collision is at "count", 
    // I use a global variable "newCount" to assign the position of the struct. 
    int newCount = count;
    for (int k = 0; k < count; k++)
    {
        for (int l = 1; (k + l) < count; l++)
        {
            Particles_Collision(&particles[k], &particles[k + l], &collision[newCount], currentTime, radius);
            newCount ++;
        }
    }


    while (currentTime < endTime)
    {
        // Do the insertion sort of collision.
        insertionSort(collision, totalCollisionStruct);

        // After insertion sort, use the collision time saved in the first struct of collision and assign it to currentTime.
        currentTime = collision[0].timeCollision;
        if (endTime < currentTime)
        {
            // If the currentTime > endTime (the time users enter), then the program will update the data struct of particle at endTime.
            for (int i = 0; i < count; i++)
            {
                particleDataUpdate(&particles[i], endTime);
            }            
        }
        else
        {
            // Use the variable, object2, to determine whether the first is Particle_Wall_Collision or Particles_Collision.
            // If object2 < 0, then it is Particle_Wall_Collision.
            if (collision[0].object2 < 0)
            {
                // Update the data of the collided particle first.
                // Then update information of other particles based on new "currentTime".  
                index = collision[0].object1;
                index2 = collision[0].object2;
                particles[index].countWall = particles[index].countWall + 1;
                particleWallUpdate(&particles[index], index2, currentTime, radius, xLength, yWidth);

                for (int i = 0; i < count; i++)
                {
                    particleDataUpdate(&particles[i], currentTime);
                }
            }
            else
            {
                // Update the data of two collided particles first.
                // Then update information of other particles based on new "currentTime".
                index = collision[0].object1;
                index2 = collision[0].object2;
                particles[index].countParticleCollision = particles[index].countParticleCollision + 1;
                particles[index2].countParticleCollision = particles[index2].countParticleCollision + 1;
                particlesUpdate(&particles[index], &particles[index2], currentTime, radius, xLength, yWidth);

                for (int i = 0; i < count; i++)
                {
                    particleDataUpdate(&particles[i], currentTime);
                }
            }
            
            // After updating the data of all particles, calculate the the new Particle_Wall_Collision and Particles_Collision.
            for (int j = 0; j < count; j++)
            {
                Particle_Wall_Collision(radius, &particles[j], &collision[j], currentTime, xLength, yWidth);
            }

            newCount = count;
            for (int k = 0; k < count; k++)
            {
                for (int l = 1; (k + l) < count; l++)
                {
                    Particles_Collision(&particles[k], &particles[(k + l)], &collision[newCount], currentTime, radius);
                    newCount = newCount + 1;
                }
                
            }
        }
    }

    // Print the information of particles.
    for (int i = 0; i < count; i++)
    {
        printf("%.6f, %.6f, %d, %d\n", particles[i].particleXPosition, particles[i].particleYPosition, particles[i].countWall, particles[i].countParticleCollision);
    }
    
    // Release the memory.
    free(particles);
    free(collision);

    return 0;
}