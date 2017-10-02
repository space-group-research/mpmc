extern "C" {
#include <mc.h>
}
// TODO: Make everything C++

// Copyright 2015 Adam Hogan

class Vec
{
 public:
  double components[3];

  Vec();

  Vec(double, double, double);

  void set(double, double, double);

  double norm();

  double dot(const Vec);

  Vec operator+(const Vec &);

  Vec operator*(const double);

  Vec &operator=(const Vec &);
};

Vec::Vec()
{
    for (int i = 0; i < 3; i++)
    {
        components[i] = 0;
    }
}

Vec::Vec(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

void Vec::set(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

double Vec::dot(const Vec other)
{
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        sum += components[i] * other.components[i];
    }
    return sum;
}

double Vec::norm()
{
    return sqrt(dot(*this));
}

Vec Vec::operator+(const Vec &right)
{
    Vec result(0, 0, 0);
    for (int i = 0; i < 3; i++)
    {
        result.components[i] = components[i] + right.components[i];
    }
    return result;
}

Vec Vec::operator*(const double x)
{
    Vec result(0, 0, 0);
    for (int i = 0; i < 3; i++)
    {
        result.components[i] = components[i] * x;
    }
    return result;
}

Vec operator*(const double x, const Vec right)
{
    Vec result(0, 0, 0);
    for (int i = 0; i < 3; i++)
    {
        result.components[i] = right.components[i] * x;
    }
    return result;
}

Vec &Vec::operator=(const Vec &right)
{
    if (&right != this)
    {  // Avoids self assignment
        components[0] = right.components[0];
        components[1] = right.components[1];
        components[2] = right.components[2];
    }
    return *this; // allows chaining, i.e. x = y = z;
}

double axilrod_teller(system_t *system)
{
    double potential = 0.0, c9, rij, rik, rjk, cos_part;
    molecule_t *molecule1, *molecule2, *molecule3;
    atom_t *atom1, *atom2, *atom3;
    pair_t temp;
    Vec ij, ik, jk, a, b;

    for (molecule1 = system->molecules; molecule1; molecule1 = molecule1->next)
    {
        for (molecule2 = system->molecules; molecule2; molecule2 = molecule2->next)
        {
            for (molecule3 = system->molecules; molecule3; molecule3 = molecule3->next)
            {
                int number_of_unique_molecules = 1;
                if (molecule2 != molecule1)
                {
                    number_of_unique_molecules++;
                }
                if (molecule3 != molecule1 && molecule3 != molecule2)
                {
                    number_of_unique_molecules++;
                }
                if (number_of_unique_molecules > 1)
                {
                    for (atom1 = molecule1->atoms; atom1; atom1 = atom1->next)
                    {
                        for (atom2 = molecule2->atoms; atom2; atom2 = atom2->next)
                        {
                            for (atom3 = molecule3->atoms; atom3; atom3 = atom3->next)
                            {
                                if (atom1 != atom2 && atom1 != atom3 && atom2 != atom3)
                                {

                                    double atom1_c9, atom2_c9, atom3_c9;
                                    atom1_c9 = atom1->c9;
                                    atom2_c9 = atom2->c9;
                                    atom3_c9 = atom3->c9;

                                    // Axilrod-Teller parameters in http://arxiv.org/pdf/1201.1532.pdf and http://dx.doi.org/10.1063/1.440310
                                    // Midzuno-Kihara approximation for c9 http://dx.doi.org/10.1143/JPSJ.11.1045
                                    if (system->midzuno_kihara_approx)
                                    {
                                        atom1_c9 = 3.0 / 4.0 * atom1->polarizability * 6.7483345 * atom1->c6;
                                        atom2_c9 = 3.0 / 4.0 * atom2->polarizability * 6.7483345 * atom2->c6;
                                        atom3_c9 = 3.0 / 4.0 * atom3->polarizability * 6.7483345 * atom3->c6;
                                    }

                                    // Mixing rule, ep. 20 of http://dx.doi.org/10.1063/1.440310
                                    c9 = pow(pow(atom1->polarizability * 6.7483345, 3) * pow(atom2->polarizability * 6.7483345, 3) * pow(atom3->polarizability * 6.7483345, 3), 1.0 / 3.0) *
                                            3.0 / (1.0 / (atom1_c9 / pow(atom1->polarizability * 6.7483345, 3)) + 1.0 / (atom2_c9 / pow(atom2->polarizability * 6.7483345, 3)) + 1.0 / (atom3_c9 / pow(atom3->polarizability * 6.7483345, 3)));

                                    if (atom1->polarizability == 0.0 || atom2->polarizability == 0.0 || atom3->polarizability == 0.0)
                                        c9 = 0.0; // avoid division by zero

                                    c9 *= 0.0032539449 / (3.166811429 * 0.000001); // convert H*Bohr^9 to K*Angstrom^9

                                    // reset fake pair ptr
                                    temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
                                    // get minimum image distance
                                    minimum_image(system, atom1, atom2, &temp);
                                    rij = temp.rimg;
                                    ij.set(temp.dimg[0], temp.dimg[1], temp.dimg[2]);

                                    // reset fake pair ptr
                                    temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
                                    // get minimum image distance
                                    minimum_image(system, atom1, atom3, &temp);
                                    rik = temp.rimg;
                                    ik.set(temp.dimg[0], temp.dimg[1], temp.dimg[2]);

                                    // reset fake pair ptr
                                    temp.d_prev[0] = temp.d_prev[1] = temp.d_prev[2] = -999999999999.;
                                    // get minimum image distance
                                    minimum_image(system, atom2, atom3, &temp);
                                    rjk = temp.rimg;
                                    jk.set(temp.dimg[0], temp.dimg[1], temp.dimg[2]);

                                    cos_part = 3;

                                    a = -1.0 * ij;
                                    b = -1.0 * ik;

                                    cos_part *= a.dot(b) / (a.norm() * b.norm());

                                    a = ij;
                                    b = -1.0 * jk;

                                    cos_part *= a.dot(b) / (a.norm() * b.norm());

                                    a = ik;
                                    b = jk;

                                    cos_part *= a.dot(b) / (a.norm() * b.norm());

                                    potential += c9 * ((1.0 + cos_part) / pow(rij * rik * rjk, 3));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // We're counting each pair 6 times
    potential = potential / 6;
    return potential;
}
