#include "quaternion.h"
#include "math.h"

/* 

©2013 Adam Hogan
Space Research Group
Department of Chemistry
University of South Florida

see http://www.cprogramming.com/tutorial/3d/quaternions.html
and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

*/

// Construct quaternion from components
void quaternion_construct_xyzw(struct quaternion *Quaternion, double x, double y, double z, double w) {
    Quaternion->x = x;
    Quaternion->y = y;
    Quaternion->z = z;
    Quaternion->w = w;
}

// Construct quaternion from an axis and an angle
// Normalizes the axis vector
void quaternion_construct_axis_angle_radian(struct quaternion *Quaternion, double x, double y, double z, double angle) {
    double magnitude = sqrt(x * x + y * y + z * z);
    if (magnitude == 0.0)  // edge case, if the axis to rotate around doesn't exist just return a quaternion with no rotation
    {
        quaternion_construct_xyzw(Quaternion, 0., 0., 0., 1.);
        return;
    }
    x = x / magnitude;
    y = y / magnitude;
    z = z / magnitude;
    double sinAngle = sin(angle / 2);
    Quaternion->x = x * sinAngle;
    Quaternion->y = y * sinAngle;
    Quaternion->z = z * sinAngle;
    Quaternion->w = cos(angle / 2);
}

// Construct quaternion from an axis and an angle (angle in degrees)
// Normalizes the axis vector
void quaternion_construct_axis_angle_degree(struct quaternion *Quaternion, double x, double y, double z, double angle) {
    angle /= 57.2957795;
    double magnitude = sqrt(x * x + y * y + z * z);
    if (magnitude == 0.0)  // edge case, if the axis to rotate around doesn't exist just return a quaternion with no rotation
    {
        quaternion_construct_xyzw(Quaternion, 0., 0., 0., 1.);
        return;
    }
    x = x / magnitude;
    y = y / magnitude;
    z = z / magnitude;
    double sinAngle = sin(angle / 2);
    Quaternion->x = x * sinAngle;
    Quaternion->y = y * sinAngle;
    Quaternion->z = z * sinAngle;
    Quaternion->w = cos(angle / 2);
}

// Normalize quaternion
void quaternion_normalize(struct quaternion *Quaternion) {
    double magnitude = sqrt(Quaternion->x * Quaternion->x + Quaternion->y * Quaternion->y + Quaternion->z * Quaternion->z + Quaternion->w * Quaternion->w);
    Quaternion->x = Quaternion->x / magnitude;
    Quaternion->y = Quaternion->y / magnitude;
    Quaternion->z = Quaternion->z / magnitude;
    Quaternion->w = Quaternion->w / magnitude;
}

// QuaternionStore = Q1 * Q2
// Order matters!
void quaternion_multiplication(struct quaternion *Q1, struct quaternion *Q2, struct quaternion *QuaternionStore) {
    double w = Q1->w * Q2->w - Q1->x * Q2->x - Q1->y * Q2->y - Q1->z * Q2->z;
    double x = Q1->w * Q2->x + Q1->x * Q2->w + Q1->y * Q2->z - Q1->z * Q2->y;
    double y = Q1->w * Q2->y - Q1->x * Q2->z + Q1->y * Q2->w + Q1->z * Q2->x;
    double z = Q1->w * Q2->z + Q1->x * Q2->y - Q1->y * Q2->x + Q1->z * Q2->w;
    QuaternionStore->w = w;
    QuaternionStore->x = x;
    QuaternionStore->y = y;
    QuaternionStore->z = z;
}

// A conjugate quaternion performs the opposite rotation
void quaternion_conjugate(struct quaternion *Quaternion, struct quaternion *QuaternionStore) {
    QuaternionStore->x = -Quaternion->x;
    QuaternionStore->y = -Quaternion->y;
    QuaternionStore->z = -Quaternion->z;
    QuaternionStore->w = Quaternion->w;
}
