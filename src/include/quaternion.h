/* 

©2013 Adam Hogan
Space Research Group
Department of Chemistry
University of South Florida

see http://www.cprogramming.com/tutorial/3d/quaternions.html
and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

*/

struct quaternion {
    double x;
    double y;
    double z;
    double w;
};

void quaternion_construct_xyzw(struct quaternion *Quaternion, double x, double y, double z, double w);
void quaternion_construct_axis_angle_radian(struct quaternion *Quaternion, double x, double y, double z, double angle);
void quaternion_construct_axis_angle_degree(struct quaternion *Quaternion, double x, double y, double z, double angle);
void quaternion_normalize(struct quaternion *Quaternion);
void quaternion_multiplication(struct quaternion *Q1, struct quaternion *Q2, struct quaternion *QuaternionStore);
void quaternion_conjugate(struct quaternion *Quaternion, struct quaternion *QuaternionStore);
